"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
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
    
    # 1. Check for the myo-inositol core with correct stereochemistry using SMARTS
    myo_inositol_pattern = Chem.MolFromSmarts('[C@@H]1([C@H]([C@@H]([C@H]([C@@H]([C@H]1O)O)O)O)O)O')
    if not mol.HasSubstructMatch(myo_inositol_pattern):
         return False, "No myo-inositol core found"

    # 2. Check for at least one phosphate group attached to the inositol
    phosphate_pattern = Chem.MolFromSmarts('[OX1,OX2][P](=[OX1])([OX1,OX2])([OX1,OX2])')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)

    if not phosphate_matches:
        return False, "No phosphate group found"
    
    return True, "Contains myo-inositol core with at least one phosphate group"