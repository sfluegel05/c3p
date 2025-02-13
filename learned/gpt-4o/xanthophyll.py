"""
Classifies: CHEBI:27325 xanthophyll
"""
from rdkit import Chem

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is an oxygenated carotenoid, characterized by a long chain
    of conjugated double bonds and one or more oxygen-containing functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for long conjugated double bonds chain (e.g., C=C-C=C-C=C)
    conjugated_pattern = Chem.MolFromSmarts("C=C.C=C.C=C.C=C")
    match_conjugated = mol.HasSubstructMatch(conjugated_pattern)

    # Check for oxygen atoms presence (common in hydroxyl, carbonyl, or ether groups in xanthophylls)
    # Basic check for presence of oxygen-containing groups
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if not match_conjugated:
        return False, "Lacks the extensive conjugated system typical for carotenoids"
    if oxygen_count < 1:
        return False, "No oxygen atoms found, required for xanthophyll classification"

    # Check for a reasonable number of double bonds term (xanthophyll typically > 9)
    double_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)
    if double_bond_count < 9:
        return False, f"Insufficient number of conjugated double bonds; found {double_bond_count}"
    
    return True, "Chemical structure fits the typical xanthophyll characteristics of conjugated double bonds and oxygen functionality"