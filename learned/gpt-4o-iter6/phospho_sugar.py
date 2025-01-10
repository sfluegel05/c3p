"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a structure with a sugar-like hydroxy pattern
    and one or more phosphate groups esterified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a phospho sugar, False otherwise.
        str: Reason for classification
    """

    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Patterns for cyclic sugars (furanose and pyranose rings)
    furanose_pattern = Chem.MolFromSmarts("C1OC[C@H]([C@@H](O)C1)O")
    pyranose_pattern = Chem.MolFromSmarts("C1O[C@H]([C@@H](CO)O)C[C@H]1O")
    
    # Patterns for linear sugars like glyceraldehyde
    linear_sugar_pattern = Chem.MolFromSmarts("CO[C@H](O)C=O") 

    # Pattern for phosphate groups (includes possibilities of di- and tri-phosphate linkages)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")

    # Check for presence of sugar structures
    has_sugar_structure = (
        mol.HasSubstructMatch(furanose_pattern) or
        mol.HasSubstructMatch(pyranose_pattern) or
        mol.HasSubstructMatch(linear_sugar_pattern)
    )
    
    if not has_sugar_structure:
        return False, "No sugar-like structure found"

    # Check for presence of phosphate group
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if not has_phosphate:
        return False, "No phosphate group esterified found"

    # If both structures are present, it's a phospho sugar
    return True, "Contains a sugar-like structure with phosphate esterification"

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