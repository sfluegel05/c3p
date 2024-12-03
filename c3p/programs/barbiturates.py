"""
Classifies: CHEBI:22693 barbiturates
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_barbiturates(smiles: str):
    """
    Determines if a molecule is a barbiturate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a barbiturate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Check for the presence of pyrimidine-2,4,6(1H,3H,5H)-trione core structure
    barbiturate_pattern = Chem.MolFromSmarts('O=C1NC(=O)NC(=O)C1')
    if not mol.HasSubstructMatch(barbiturate_pattern):
        return False, "Pyrimidine-2,4,6(1H,3H,5H)-trione core structure not found"

    return True, "Molecule is a barbiturate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22693',
                          'name': 'barbiturates',
                          'definition': 'Members of the class of pyrimidones '
                                        'consisting of '
                                        'pyrimidine-2,4,6(1H,3H,5H)-trione '
                                        '(barbituric acid) and its '
                                        'derivatives. Largest group of the '
                                        'synthetic sedative/hypnotics, sharing '
                                        'a characteristic six-membered ring '
                                        'structure.',
                          'parents': ['CHEBI:38337']},
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
    'num_true_positives': 9,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9,
    'f1': 0.9473684210526316,
    'accuracy': None}