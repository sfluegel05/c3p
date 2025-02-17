"""
Classifies: CHEBI:47908 alkanethiol
"""
"""
Classifies: alkanethiol
Definition: An alkanethiol is a compound in which a sulfanyl group (-SH) is attached to an alkyl group.
Examples include 3-Methyl-2-butanethiol, propane-2-thiol, ethanethiol, etc.
"""

from rdkit import Chem

def is_alkanethiol(smiles: str):
    """
    Determines if a molecule is an alkanethiol based on its SMILES string.
    An alkanethiol is defined as a compound with a thiol (–SH) group attached to an alkyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an alkanethiol, False otherwise
        str: Explanation for the classification decision
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update property cache and compute implicit hydrogens
    mol.UpdatePropertyCache()
    
    # Flag to indicate if any valid alkanethiol group is found.
    valid_thiol_found = False
    invalid_thiol_found = False

    # Iterate over atoms looking for sulfur atoms with thiol properties
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # sulfur atom
            # Check for at least one hydrogen attached. In RDKit, hydrogens may be implicit.
            # GetTotalNumHs() returns the total number of hydrogens (implicit+explicit)
            if atom.GetTotalNumHs() < 1:
                continue  # not a thiol since no hydrogen on sulfur

            # Get heavy-atom neighbors (neighbors that are not hydrogen, as implicit H are not present in GetNeighbors())
            neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
            # For a –SH group, we typically expect the S to be attached to exactly one heavy atom.
            if len(neighbors) != 1:
                # This sulfur might be in another context (e.g. bridging sulfur or disulfide) so we mark it 
                invalid_thiol_found = True
                continue

            alkyl_candidate = neighbors[0]
            # The neighbor must be a carbon (atomic number 6)
            if alkyl_candidate.GetAtomicNum() != 6:
                invalid_thiol_found = True
                continue

            # Now check that the carbon is saturated (sp3 hybridized) and not aromatic.
            if alkyl_candidate.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or alkyl_candidate.GetIsAromatic():
                invalid_thiol_found = True
                continue

            # If we passed all these tests then we have found a thiol group attached to an alkyl chain.
            valid_thiol_found = True
            # We can stop early since one valid group is sufficient.
            break

    if valid_thiol_found:
        return True, "Contains at least one –SH group attached to a saturated alkyl (sp3 carbon) chain."
    elif invalid_thiol_found:
        return False, "Thiolate group found but it is not attached to a saturated alkyl (sp3 carbon) group."
    else:
        return False, "No thiol (–SH) group found in the molecule."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "SC(C(C)C)C",   # 3-Methyl-2-butanethiol
        "CCS",          # ethanethiol
        "SCCCC",        # butanethiol
        "S\C=C(\CC)/C"  # 2-Methyl-1-propenethiol (should return False because attached carbon is sp2)
    ]
    
    for smi in test_smiles:
        is_valid, reason = is_alkanethiol(smi)
        print(f"SMILES: {smi:20s} -> {is_valid} ({reason})")