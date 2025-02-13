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
        bool: True if the molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Redefining a pattern for sugar-like hydroxy structures
    # This accounts for both furanose and pyranose ring sugars
    sugar_ring_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C1O")  # Typical hexose as an example
    # Linear sugar parts (e.g., glyceraldehyde phosphate)
    sugar_linear_pattern = Chem.MolFromSmarts("OCC(O)CO")

    # Refine patterns for phosphate groups (covering monophosphate, diphosphate, etc.)
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)([O])[O]")  # Phosphate group

    # Check for presence of sugar structure
    has_sugar_structure = mol.HasSubstructMatch(sugar_ring_pattern) or mol.HasSubstructMatch(sugar_linear_pattern)
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