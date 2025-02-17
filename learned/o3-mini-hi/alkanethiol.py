"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: An alkanethiol is a compound in which a sulfanyl (-SH) group is attached to an alkyl group.
This version excludes molecules that contain peptide (amide) bonds and other functionalities (e.g. carboxyl, sulfonic, phosphate)
that would indicate the thiol is embedded in a multifunctional (often biomolecular) environment.
It also accepts cases where the –SH is attached to a non‐aromatic carbon even if that carbon is sp2 hybridized.
"""

from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as a compound with a thiol (–SH) group attached 
    to an alkyl group (non-aromatic carbon); the definition here also rejects molecules 
    with extra functionalities (e.g. amide, carboxyl, sulfonic or phosphate groups) that suggest
    the thiol is part of a biomolecule or a more complex system.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as an alkanethiol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol.UpdatePropertyCache()
    
    # Exclude molecules with amide bonds (common in peptides)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_smarts):
        return False, "Molecule contains amide bonds (suggesting a peptide backbone)"
    
    # Exclude molecules containing carboxyl groups together with amine functions (e.g. amino acids)
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
    
    # Iterate over all sulfur atoms and check for –SH groups.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # Not a sulfur atom
        
        # Check if this sulfur atom has at least one hydrogen attached.
        # (This works even if H's are implicit.)
        if atom.GetTotalNumHs() < 1:
            continue  # Not a thiol if no hydrogen is attached
        
        # Get heavy-atom neighbors (ignore hydrogens)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        
        # For a genuine thiol group, we expect exactly one heavy neighbor.
        if len(heavy_neighbors) != 1:
            invalid_thiol_found = True
            continue
        
        neighbor = heavy_neighbors[0]
        # The neighbor should be a carbon (atomic number 6) to be an alkyl linkage.
        if neighbor.GetAtomicNum() != 6:
            invalid_thiol_found = True
            continue
        
        # We now relax the hybridization check: instead of requiring sp3 we only disallow aromatic atoms.
        if neighbor.GetIsAromatic():
            invalid_thiol_found = True
            continue

        # Additional check: if the neighbor is directly attached to a carbonyl carbon (C=O) then it may be part of a carbonyl functional group.
        for nbr in neighbor.GetNeighbors():
            if nbr.GetIdx() == atom.GetIdx():
                continue
            # Check if neighbor has a double-bonded oxygen.
            for bond in mol.GetBondBetweenAtoms(neighbor.GetIdx(), nbr.GetIdx()):
                if bond.GetBondTypeAsDouble() == 2 and nbr.GetAtomicNum() == 8:
                    invalid_thiol_found = True
                    break
            if invalid_thiol_found:
                break
        if invalid_thiol_found:
            continue

        # If we got here, we found a candidate –SH attached to a non‐aromatic carbon (even if sp2)
        valid_thiol_found = True
        # If at least one valid thiol is found, we consider the molecule an alkanethiol.
        break

    if valid_thiol_found:
        return True, "Contains at least one –SH group attached to a non-aromatic alkyl carbon."
    elif invalid_thiol_found:
        return False, "Thiolate group present but attached carbon is either aromatic, conjugated with a carbonyl, or otherwise inconsistent with a simple alkyl group."
    else:
        return False, "No –SH (thiol) group found in the molecule."

# Example usage when running this script directly.
if __name__ == "__main__":
    test_smiles = [
        "SC(C(C)C)C",            # 3-Methyl-2-butanethiol (true positive)
        "CC(C)S",                # propane-2-thiol (true positive)
        "SC(C(C)C)CO",           # xi-2-Mercapto-3-methyl-1-butanol (true positive)
        "SC(CCC)C",              # 2-Pentanethiol (true positive)
        "SC(C(CO)C)CC",          # 3-Mercapto-2-methylpentanol (true positive)
        "SCCCCCCCCCS",           # 1,9-Nonanedithiol (true positive)
        "CCS",                   # ethanethiol (true positive)
        "Cl.SCCN",               # Cysteamine hydrochloride (true positive)
        "S\C=C(\CC)/C",          # 2-Methyl-1-propenethiol (should be classified as alkanethiol now)
        "NC(CCS)C(O)=O",         # homocysteine (false positive, rejected)
    ]
    
    for smi in test_smiles:
        result, explanation = is_alkanethiol(smi)
        print(f"SMILES: {smi:50s} -> {result} ({explanation})")