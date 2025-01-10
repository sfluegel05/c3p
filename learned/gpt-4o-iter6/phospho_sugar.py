"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a polyhydroxy backbone typical of monosaccharides with at least one phosphate group esterified.

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

    # Define a broad pattern for polyhydroxy structures typical of monosaccharides
    polyol_pattern = Chem.MolFromSmarts("[C,O]([C,O])([O][H])[C,O]([C,O])([O][H])")

    # Define patterns for phosphate groups (both mono- and diesterified forms)
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)([O,H])([O,H])")

    # Check for presence of polyhydroxy structure
    has_polyhydroxy = mol.HasSubstructMatch(polyol_pattern)
    if not has_polyhydroxy:
        return False, "No polyhydroxy structure typical of monosaccharides found"

    # Check for presence of phosphate group
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if not has_phosphate:
        return False, "No phosphate group esterified found"

    # If both structures are present, it's a phospho sugar
    return True, "Contains a polyhydroxy backbone with phosphate esterification"

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