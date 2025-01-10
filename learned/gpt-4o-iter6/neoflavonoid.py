"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define benzopyran ring system SMARTS and 4-aryl substitution
    benzopyran_pattern = Chem.MolFromSmarts("c1ccccc1O2CC=C2")  # Simplified 1-benzopyran
    aryl_at_pos_4_pattern = Chem.MolFromSmarts("c1ccccc1")  # Aryl group

    # Check for benzopyran backbone structure
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 1-benzopyran structure found"

    # Find the position of the benzopyran in the molecule
    benzopyran_matches = mol.GetSubstructMatch(benzopyran_pattern)
    
    if not benzopyran_matches:
        return False, "1-benzopyran structure not properly matched in the molecule"

    # Assume fourth atom in benzopyran_matches is the position of interest
    pos_4_atom_idx = benzopyran_matches[3] # 4th position based on 1-benzopyran
    
    # Check if an aryl group is attached here
    for neighbor in mol.GetAtomWithIdx(pos_4_atom_idx).GetNeighbors():
        if mol.HasSubstructMatch(aryl_at_pos_4_pattern, atoms=[neighbor.GetIdx()]):
            return True, "Aryl substituent found at position 4"
    
    return False, "No aryl substitution at position 4 in the 1-benzopyran"