"""
Classifies: CHEBI:26401 purines
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_purines(smiles: str):
    """
    Determines if a molecule is a purine (purine or substituted purine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a purine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for purine
    purine_smarts = "c1ncnc2ncnc12"
    purine_pattern = Chem.MolFromSmarts(purine_smarts)

    if mol.HasSubstructMatch(purine_pattern):
        return True, "Molecule contains the purine core structure"
    else:
        return False, "Molecule does not contain the purine core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26401',
                          'name': 'purines',
                          'definition': 'A class of imidazopyrimidines that '
                                        'consists of purine and its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:35875']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 165,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9939759036144579,
    'f1': 0.9969788519637462,
    'accuracy': None}