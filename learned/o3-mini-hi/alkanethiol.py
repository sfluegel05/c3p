"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: An alkanethiol is a compound in which a sulfanyl (-SH) group is attached to an alkyl group.
This implementation excludes molecules that contain peptide (amide) bonds or additional functionalities
(e.g., carboxyl, sulfonic, phosphate) that might indicate the –SH is embedded in a more complex biomolecular framework.
"""

from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as a compound with a thiol (-SH) group attached 
    to an alkyl group. Extra functionalities (e.g., amide, carboxyl, sulfonic or phosphate groups)
    cause the molecule to be rejected as a simple alkanethiol.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an alkanethiol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol.UpdatePropertyCache()
    
    # Exclude molecules with amide bonds (common in peptides)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_smarts):
        return False, "Molecule contains amide bonds (suggests peptide backbone)"
    
    # Exclude molecules that contain both carboxyl and amine groups (e.g., amino acids)
    carboxyl_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1,H0,-1]")
    amine_smarts = Chem.MolFromSmarts("[NX3;H2,H1]")
    if mol.HasSubstructMatch(carboxyl_smarts) and mol.HasSubstructMatch(amine_smarts):
        return False, "Molecule contains both amine and carboxyl groups (likely a biofunctional molecule)"
    
    # Exclude molecules with sulfonic acid or phosphoryl groups
    sulfonic_smarts = Chem.MolFromSmarts("[SX4](=O)(=O)[OX1,H0,-1]")
    phosphate_smarts = Chem.MolFromSmarts("P(=O)(O)O")
    if mol.HasSubstructMatch(sulfonic_smarts):
        return False, "Molecule contains sulfonic acid functionalities not expected for a simple alkanethiol"
    if mol.HasSubstructMatch(phosphate_smarts):
        return False, "Molecule contains phosphate groups not expected for a simple alkanethiol"
    
    valid_thiol_found = False
    invalid_thiol_found = False
    
    # Loop over all atoms looking for sulfur atoms that might be part of a thiol group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # Not a sulfur atom
        
        # Check if the sulfur atom has at least one hydrogen attached (explicitly or implicitly).
        if atom.GetTotalNumHs() < 1:
            continue  # Not a thiol if no hydrogen is attached
        
        # Get heavy-atom (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        
        # For a genuine thiol group, expect exactly one heavy neighbor.
        if len(heavy_neighbors) != 1:
            invalid_thiol_found = True
            continue
        
        neighbor = heavy_neighbors[0]
        # The neighbor should be carbon (atomic number 6) to be considered an alkyl linkage.
        if neighbor.GetAtomicNum() != 6:
            invalid_thiol_found = True
            continue
        
        # Exclude cases where the neighbor is aromatic.
        if neighbor.GetIsAromatic():
            invalid_thiol_found = True
            continue

        # Check that the neighbor is not directly attached to a carbonyl group (C=O).
        for nbr in neighbor.GetNeighbors():
            # Skip back edge to the sulfur atom.
            if nbr.GetIdx() == atom.GetIdx():
                continue
            bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx())
            if bond is not None:
                if bond.GetBondTypeAsDouble() == 2 and nbr.GetAtomicNum() == 8:
                    invalid_thiol_found = True
                    break
        if invalid_thiol_found:
            continue

        # Found a valid –SH group attached to an alkyl (non-aromatic) carbon.
        valid_thiol_found = True
        break

    if valid_thiol_found:
        return True, "Contains at least one –SH group attached to a non-aromatic alkyl carbon."
    elif invalid_thiol_found:
        return False, "Thiolate group present but attached carbon is either aromatic, directly conjugated with a carbonyl, or otherwise inconsistent with a simple alkyl group."
    else:
        return False, "No valid –SH (thiol) group found in the molecule."

# Example usage when running this script directly.
if __name__ == "__main__":
    test_smiles = [
        "SC(C(C)C)C",            # 3-Methyl-2-butanethiol
        "CC(C)S",                # propane-2-thiol
        "SC(C(C)C)CO",           # xi-2-Mercapto-3-methyl-1-butanol
        "SC(CCC)C",              # 2-Pentanethiol
        "SC(C(CO)C)CC",          # 3-Mercapto-2-methylpentanol
        "SC(CC)CS",              # 1,2-Butanedithiol
        "SCCCCCCCCCS",           # 1,9-Nonanedithiol
        "CCS",                   # ethanethiol
        "Cl.SCCN",               # Cysteamine hydrochloride
        "S\\C=C(\\CC)/C",        # 2-Methyl-1-propenethiol
        "CC(=O)CC(C)(C)S",       # 4-mercapto-4-methylpentan-2-one
        "NCCCNCCS",              # WR-1065
    ]
    
    for smi in test_smiles:
        result, explanation = is_alkanethiol(smi)
        print(f"SMILES: {smi:30s} -> {result} ({explanation})")