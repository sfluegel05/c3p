"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is an aldehyde where the carbonyl group is directly attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized check for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("[c]1:[c]:[c]:[c]:[c]:[c]1")
    aromatic_pattern2 = Chem.MolFromSmarts("[c]1:[c]:[c]:[c]:[c]:[n]1")
    aromatic_pattern3 = Chem.MolFromSmarts("[c]1:[c]:[c]:[c]:[o]:[c]1")
    aromatic_pattern4 = Chem.MolFromSmarts("[c]1:[c]:[c]:[c]:[n]:[n]1")
    aromatic_pattern5 = Chem.MolFromSmarts("[c]1:[c]:[c]:[c]:[n]:[c]1")
    aromatic_pattern6 = Chem.MolFromSmarts("[c]1:[c]:[c]:[c]:[n]:[o]1")
    aromatic_pattern7 = Chem.MolFromSmarts("[c]1:[c]:[c]:[s]:[c]:[c]1")
    has_aromatic = mol.HasSubstructMatch(aromatic_pattern) or \
                   mol.HasSubstructMatch(aromatic_pattern2) or \
                   mol.HasSubstructMatch(aromatic_pattern3) or \
                   mol.HasSubstructMatch(aromatic_pattern4) or \
                   mol.HasSubstructMatch(aromatic_pattern5) or\
                   mol.HasSubstructMatch(aromatic_pattern6) or \
                   mol.HasSubstructMatch(aromatic_pattern7)


    if not has_aromatic:
      return False, "No aromatic ring found"

    # Check for aldehyde group directly attached to aromatic carbon
    aldehyde_pattern = Chem.MolFromSmarts("[cH0,cH1][CX3](=[OX1])[#6H0]")

    if not mol.HasSubstructMatch(aldehyde_pattern):
       return False, "No aldehyde group directly attached to aromatic ring found"

    return True, "Contains an aromatic ring with a directly attached aldehyde group."