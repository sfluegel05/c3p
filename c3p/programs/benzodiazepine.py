"""
Classifies: CHEBI:22720 benzodiazepine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_benzodiazepine(smiles: str):
    """
    Determines if a molecule is a benzodiazepine (benzene ring fused to a diazepine ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzodiazepine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for benzodiazepine
    benzodiazepine_smarts = "c1ccccc1C2=NCCN=C2"  # Simplified SMARTS for benzodiazepine core structure
    pattern = Chem.MolFromSmarts(benzodiazepine_smarts)

    if mol.HasSubstructMatch(pattern):
        return True, "Molecule contains the benzodiazepine core structure"
    else:
        return False, "Molecule does not contain the benzodiazepine core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22720',
                          'name': 'benzodiazepine',
                          'definition': 'A group of heterocyclic compounds '
                                        'with a core structure containing a '
                                        'benzene ring fused to a diazepine '
                                        'ring.',
                          'parents': ['CHEBI:38101', 'CHEBI:38166']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 17,
    'num_false_negatives': 17,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}