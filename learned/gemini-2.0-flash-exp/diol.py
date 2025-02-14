"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a molecule with at least two hydroxyl (-OH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None

    # Count hydroxyl groups
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8: # Check for oxygen atoms
            for neighbor in atom.GetNeighbors(): # Check oxygen neighbors
                if neighbor.GetAtomicNum() == 1: # Check for hydrogen neighbors
                    hydroxyl_count += 1 # Increment if oxygen-hydrogen bond found
                    break
    
    if hydroxyl_count >= 2:
        return True, f"Molecule contains at least two hydroxyl groups. It has {hydroxyl_count}."
    else:
       return False, f"Molecule has fewer than two hydroxyl groups, it has {hydroxyl_count}."