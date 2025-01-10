"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of hydroxyl groups (-OH)
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:  # Hydrogen
                    hydroxyl_count += 1

    # Check if the molecule has exactly two hydroxyl groups
    if hydroxyl_count == 2:
        return True, "Contains exactly two hydroxyl groups"
    else:
        return False, f"Contains {hydroxyl_count} hydroxyl groups, need exactly 2"