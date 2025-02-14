"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is a myo-inositol with at least one phosphate group directly attached to the ring carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the myo-inositol core with correct stereochemistry
    myo_inositol_pattern = Chem.MolFromSmarts('[C@@H]1([C@H]([C@@H]([C@H]([C@@H]([C@H]1O)O)O)O)O)O')
    if not mol.HasSubstructMatch(myo_inositol_pattern):
         return False, "No myo-inositol core found"

    # 2. Check that all the 6 carbons in ring are part of the inositol core
    c_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(f'[C@@H]1([C@H]([C@@H]([C@H]([C@@H]([C@H]1O)O)O)O)O)O')):
                 c_count +=1
    if c_count != 6:
        return False, "Not all carbons in the six-membered ring are part of the inositol"

    # 3. Check for at least one phosphate group directly attached to a carbon of the inositol ring using SMARTS.
    # Use a more specific SMARTS pattern to check for phosphate attached directly to the carbons
    phosphate_pattern = Chem.MolFromSmarts("[C;R][O][P](=[O])(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
         return False, "No phosphate group directly attached to the myo-inositol ring found"

    # 4. Exclude molecules with lipid or sugar attachments
    # Check for glycerol (lipid) attachments 
    glycerol_pattern = Chem.MolFromSmarts("C(O)(C(O)C(O))")
    if mol.HasSubstructMatch(glycerol_pattern):
        return False, "Molecule has a glycerol lipid component, not a free myo-inositol phosphate"
    
    # Check for sugar attachments (excluding inositol itself)
    sugar_pattern = Chem.MolFromSmarts('OC[C,O][C,O][C,O][C,O][C,O][C,O]')
    if mol.HasSubstructMatch(sugar_pattern):
        sugar_matches = mol.GetSubstructMatches(sugar_pattern)
        is_inositol = False
        for match in sugar_matches:
            if set([atom.GetIdx() for atom in mol.GetAtoms()]) == set(match):
               is_inositol = True
        if not is_inositol:
            return False, "Molecule has a sugar component, not a free myo-inositol phosphate"



    
    return True, "Contains myo-inositol core with at least one phosphate group directly attached to the ring carbons"