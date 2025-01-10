"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion is a derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen have been replaced with univalent (usually organyl) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define quaternary ammonium ion SMARTS pattern
    # Nitrogen atom with formal charge +1 and degree 4 (four single bonds)
    quat_ammonium_pattern = Chem.MolFromSmarts("[N+;D4]")
    if mol.HasSubstructMatch(quat_ammonium_pattern):
        return True, "Contains quaternary ammonium ion"
    else:
        return False, "No quaternary ammonium ion found"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'quaternary ammonium ion',
        'definition': 'A derivative of ammonium, NH4(+), in which all four of the hydrogens bonded to nitrogen have been replaced with univalent (usually organyl) groups.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}