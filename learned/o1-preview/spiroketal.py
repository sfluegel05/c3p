"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:35567 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for spiroketal
    # [C;R2;X4] - spiro carbon atom in two rings with valence 4
    # Attached to two oxygens in rings
    spiroketal_smarts = "[C;R2;X4](-[O;R])(-[O;R])(-[*;R])(-[*;R])"
    spiroketal_pattern = Chem.MolFromSmarts(spiroketal_smarts)

    # Find matches in the molecule
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    if matches:
        spiro_idx = matches[0][0]  # Get the index of the spiro carbon atom
        return True, f"Spiroketal detected at atom index {spiro_idx}"
    else:
        return False, "No spiroketal center found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35567',
                              'name': 'spiroketal',
                              'definition': 'A cyclic ketal in which the ketal '
                                            'carbon is the only common atom of '
                                            'two rings.',
                              'parents': []},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}