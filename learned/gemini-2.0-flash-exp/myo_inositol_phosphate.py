"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

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

    # 2. Check for at least one phosphate group directly attached to a carbon of the inositol ring using SMARTS.
    phosphate_pattern = Chem.MolFromSmarts("[C;R][O][P](=[O])(O)(O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
         return False, "No phosphate group directly attached to the myo-inositol ring found"

    return True, "Contains myo-inositol core with at least one phosphate group directly attached to the ring carbons"