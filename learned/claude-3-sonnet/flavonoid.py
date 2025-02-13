"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:35508 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid is a compound whose skeleton is based on 1-benzopyran with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define flavonoid core pattern
    flavonoid_core = Chem.MolFromSmarts('[o;r5]1[c;r5]2[c;r6]([c;r5][c;r5][c;r6]1)[c;r6][c;r6][c;r6]2[c;r5]3[c;r6][c;r6][c;r5][c;r6][c;r6]3')
    
    # Check for flavonoid core
    match = mol.GetSubstructMatch(flavonoid_core)
    if not match:
        return False, "Missing flavonoid core skeleton"
    
    # Check for aryl substituent at position 2
    aryl_sub_pattern = Chem.MolFromSmarts('[c;r6][a;r6][a;r6][a;r6]')
    aryl_sub_match = mol.GetSubstructMatches(aryl_sub_pattern)
    for idx in aryl_sub_match:
        if mol.GetBondBetweenAtoms(idx[0], match[2]).GetBondType() == Chem.BondType.SINGLE:
            break
    else:
        return False, "No aryl substituent found at position 2"
    
    # Check for common flavonoid substituents
    sub_patterns = [
        Chem.MolFromSmarts('[OH]'),      # hydroxy
        Chem.MolFromSmarts('[OC]'),      # methoxy
        Chem.MolFromSmarts('[O;r6]'),    # glycoside
        Chem.MolFromSmarts('[C=C]'),     # prenyl/isoprenyl
        Chem.MolFromSmarts('[CX3](=O)'), # acyl
    ]
    sub_matches = [mol.GetSubstructMatches(pat) for pat in sub_patterns]
    if not any(sub_matches):
        return False, "No common flavonoid substituents found"
    
    return True, "Contains flavonoid core skeleton with aryl substituent at position 2"