"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string. 
    An aliphatic alcohol has a hydroxyl group (-OH) attached to a non-aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for aliphatic alcohol (non-aromatic C attached to OH)
    aliphatic_oh_pattern = Chem.MolFromSmarts("[C;!R][OX2H]")  # C not in a ring and OH group

    # Find all OH groups that satisfy aliphatic pattern
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            if not atom.GetIsAromatic():  # Ensuring the oxygen is not part of aromatic system
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and not neighbor.GetIsAromatic():
                        return True, "Contains aliphatic hydroxyl group(s) on a non-aromatic carbon chain"

    return False, "No aliphatic hydroxyl group found"