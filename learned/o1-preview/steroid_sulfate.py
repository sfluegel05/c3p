"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a steroid molecule where at least one hydroxy group is esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone pattern (cyclopentanoperhydrophenanthrene nucleus)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3C(C4C=C5)=C5C4C3CCC2C1')  # Adjusted SMARTS pattern for steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Define sulfate ester group pattern
    sulfate_pattern = Chem.MolFromSmarts('O[S](=O)(=O)[O]')
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    if not sulfate_matches:
        return False, "No sulfate ester group found"

    # Check that sulfate group is attached via oxygen to the steroid backbone
    steroid_matches = mol.GetSubstructMatch(steroid_pattern)
    steroid_atom_indices = set(steroid_matches)

    # For each sulfate ester group, check connection to steroid backbone
    for match in sulfate_matches:
        # The ester oxygen atom connected to sulfur
        ester_oxygen_idx = match[0]
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
        # Check neighbors of ester oxygen atom
        for neighbor in ester_oxygen_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in steroid_atom_indices:
                return True, "Contains steroid backbone with sulfate ester group attached via oxygen"

    return False, "Sulfate group not attached to steroid backbone via ester linkage"


__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'steroid sulfate',
        'definition': 'A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.',
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
    'stdout': None
}