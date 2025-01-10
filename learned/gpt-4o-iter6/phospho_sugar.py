"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for monosaccharide ring(s) with hydroxy groups
    # We assume rings of sizes 5 or 6 are required (common for sugars)
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No ring structures found, expected monosaccharide backbone"
    
    # Looking for a structure consistent with monosaccharides with hydroxy groups
    sugar_pattern = Chem.MolFromSmarts("CO[C@H]1O[C@@H]1")  # Simplified monosaccharide pattern
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No monosaccharide backbone with hydroxy groups found"
    
    # Looking for esterified phosphate group (C-O-P(=O)(O)-O)
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester groups found"
    
    return True, "Contains a monosaccharide backbone with hydroxy groups esterified with phosphoric acid"

# Example metadata for the classification
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:??',
        'name': 'phospho_sugar',
        'definition': 'Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.'
    },
    'message': None,
    'success': True
}