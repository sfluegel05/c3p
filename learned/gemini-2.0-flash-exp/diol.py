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
    o_count = 0
    h_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
           o_count += 1
        elif atom.GetAtomicNum() == 1:
           h_count +=1

    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                  hydroxyl_count+=1
    
    if hydroxyl_count >= 2 and hydroxyl_count == o_count and hydroxyl_count ==h_count:
        return True, f"Molecule contains at least two hydroxyl groups. It has {hydroxyl_count}."
    elif hydroxyl_count < 2:
       return False, f"Molecule has fewer than two hydroxyl groups, it has {hydroxyl_count}."
    else:
       return False, f"Molecule has {hydroxyl_count} O-H but they don't all match"