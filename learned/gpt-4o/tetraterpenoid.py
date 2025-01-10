"""
Classifies: CHEBI:26935 tetraterpenoid
"""
from rdkit import Chem

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    
    Tetraterpenoids are derived from a C40 backbone, involving polyene structures, possibly consisting of cyclic forms, oxygens, and allowing C-skeleton variations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check number of carbon atoms, providing flexibility for isomeric/diverse structures
    num_carbons = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if num_carbons < 30 or num_carbons > 55:  # Broader range to account for diversity
        return False, f"Unusual number of carbons for tetraterpenoid: {num_carbons} carbons"

    # Check presence of multiple double bonds characteristic of polyene structures
    num_double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() >= 2)
    if num_double_bonds < 5:  # Lower threshold for flexibility
        return False, "Insufficient presence of double bonds for a tetraterpenoid structure"
    
    # Consider characteristic cyclic elements or functional groups
    # This might typically involve looking for potential ring structures or oxy-functionals
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No significant cyclic structures suggested"
    
    # Check for isoprene-like repeating units, but more flexible in counting
    isoprene_unit = Chem.MolFromSmarts("C=C-C-C=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_unit)
    if len(isoprene_matches) < 4:  # Even more relaxed threshold
        return False, "Not sufficient isoprene-like units present"
    
    # If all checks pass
    return True, "Structure likely aligns with characteristic features of a tetraterpenoid"