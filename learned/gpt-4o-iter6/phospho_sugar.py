"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains any sugar-like structure with one or more phosphate groups esterified.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Patterns for sugar-like structures
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@@H]"),  # General sugar motifs
        Chem.MolFromSmarts("OC[C@H](O)C(O)"),  # Specific for furanose and pyranose
        Chem.MolFromSmarts("[C@H]1O[C@@H](O)[C@@H](O)[C@H]1"),  # Ring sugars
        Chem.MolFromSmarts("CO[C@H](O)C")  # Linear sugar motifs
    ]

    # Patterns for mono-, di-, and triphosphate groups
    phosphate_patterns = [
        Chem.MolFromSmarts("O=P(O)(O)O"),         # Monophosphate
        Chem.MolFromSmarts("O=P([O-])(=O)OP([O-])(=O)O"), # Diphosphate
        Chem.MolFromSmarts("O=P([O-])(=O)OP([O-])(=O)OP([O-])(=O)O")  # Triphosphate
    ]

    # Check for the presence of sugar and phosphate structures
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    has_phosphate = any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns)

    if not has_sugar:
        return False, "No recognizable sugar structure found"

    if not has_phosphate:
        return False, "No esterified phosphate group found"

    return True, "Contains a sugar structure with phosphate esterification"

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