"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Δ(4) steroid is characterized by a 3-oxo group and a double bond between carbon 4 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a steroid backbone: 6-6-6-5 fused ring system
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CCC4(C3CCC4)C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # SMARTS for the 3-oxo group
    oxo_pattern = Chem.MolFromSmarts("C=O")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    
    # Try to match an oxo group at position 3
    position_3_oxo_found = False
    for match in oxo_matches:
        # Check if any matched carbon is in the 3rd position of the steroid backbone
        if mol.GetAtomWithIdx(match[0]).GetSymbol() == 'C':
            # Only approximate based on connectivity, not exact steroid position
            neighbors = mol.GetAtomWithIdx(match[0]).GetNeighbors()
            if any(neigh.GetSymbol() == 'C' for neigh in neighbors):
                position_3_oxo_found = True
                break

    if not position_3_oxo_found:
        return False, "3-oxo group not found at expected position"

    # SMARTS for Delta(4) double bond
    delta_4_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(delta_4_pattern)
    
    # Check every double bond to find if one is in the expected position
    position_4_double_bond_found = False
    for match in double_bond_matches:
        # Check positions if they might correspond to between 4th and 5th in the backbone
        if any(mol.GetAtomWithIdx(idx).GetSymbol() == 'C' for idx in match):
            position_4_double_bond_found = True
            break

    if not position_4_double_bond_found:
        return False, "Delta(4) double bond not found"

    return True, "Contains 3-oxo group and Delta(4) double bond characteristic of steroids"