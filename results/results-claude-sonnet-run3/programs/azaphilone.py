from rdkit import Chem
from rdkit.Chem import AllChem

def is_azaphilone(smiles: str):
    """
    Determines if a molecule is an azaphilone based on the presence of a 
    6H-isochromene-6,8(7H)-dione or isoquinoline-6,8(2H,7H)-dione skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azaphilone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core structure SMARTS patterns
    # Pattern 1: Basic isochromene-dione core with specific ketone positions
    pattern1 = Chem.MolFromSmarts('[#6]1~[#6]~[#6]2~[#6](=[#8])~[#6](~[#6]2~[#8]~[#6]1)~[#6](=[#8])')
    
    # Pattern 2: Alternative isochromene-dione core with different bond arrangement
    pattern2 = Chem.MolFromSmarts('[#6]1~[#6]~[#6]2~[#6](~[#8]~[#6]1)~[#6](=[#8])~[#6]2=[#8]')
    
    # Pattern 3: Isoquinoline-dione core
    pattern3 = Chem.MolFromSmarts('[#6]1~[#6]~[#6]2~[#6](~[#6]~[#7]1)~[#6](=[#8])~[#6]2=[#8]')
    
    # Pattern 4: Alternative pattern capturing more variations
    pattern4 = Chem.MolFromSmarts('[#6]1~[#6]2~[#6](~[#6]~[#8,#7]~[#6]1)~[#6](=[#8])~[#6](~[#6]2)=[#8]')

    # Check each pattern
    if mol.HasSubstructMatch(pattern1):
        return True, "Contains isochromene-dione core structure (type 1)"
    elif mol.HasSubstructMatch(pattern2):
        return True, "Contains isochromene-dione core structure (type 2)"
    elif mol.HasSubstructMatch(pattern3):
        return True, "Contains isoquinoline-dione core structure"
    elif mol.HasSubstructMatch(pattern4):
        return True, "Contains azaphilone core structure"

    return False, "Does not contain required azaphilone core structure"
# Pr=None
# Recall=0.0