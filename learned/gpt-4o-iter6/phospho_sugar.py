"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a sugar-like structure with one or more phosphate groups esterified.
    
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
    
    # Patterns for common sugar-like structures
    cyclic_sugar_patterns = [
        "C1OC[C@H](O1)C(O)",  # furanose
        "C1O[C@H](O)C[C@H]1O",  # pyranose
    ]

    linear_sugar_patterns = [
        "CO[C@H](O)C=O"  # glyceraldehyde-like
    ]

    phosphate_patterns = [
        "O=P(O)(O)O"  # monophosphate
    ]

    # Check for presence of sugar and phosphate structures
    has_sugar = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in cyclic_sugar_patterns + linear_sugar_patterns)
    has_phosphate = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in phosphate_patterns)

    if not has_sugar:
        return False, "No recognizable sugar structure found"
    
    if not has_phosphate:
        return False, "No esterified phosphate group found"

    # Confirm sugar-phosphate linkage
    sugar_phosphate_linkage = Chem.MolFromSmarts("CO[P](=O)(O)O")
    if not mol.HasSubstructMatch(sugar_phosphate_linkage):
        return False, "No sugar-phosphate linkage found"

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