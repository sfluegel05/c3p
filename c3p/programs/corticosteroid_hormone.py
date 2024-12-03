"""
Classifies: CHEBI:36699 corticosteroid hormone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_corticosteroid_hormone(smiles: str):
    """
    Determines if a molecule is a corticosteroid hormone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a corticosteroid hormone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the steroid backbone
    steroid_core = Chem.MolFromSmarts('C1CCC2C1CCC3C2CCC4C3CC(O)CC4')
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Molecule does not have the steroid backbone"

    # Check for functional groups characteristic of corticosteroids
    functional_groups = [
        Chem.MolFromSmarts('C=O'),  # Ketone group
        Chem.MolFromSmarts('O'),    # Hydroxyl group
        Chem.MolFromSmarts('F'),    # Fluorine atom
        Chem.MolFromSmarts('Cl'),   # Chlorine atom
        Chem.MolFromSmarts('Br')    # Bromine atom
    ]

    for fg in functional_groups:
        if mol.HasSubstructMatch(fg):
            return True, "Molecule has the steroid backbone and characteristic functional groups of corticosteroids"

    return False, "Molecule does not have the characteristic functional groups of corticosteroids"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36699',
                          'name': 'corticosteroid hormone',
                          'definition': 'Any of a class of steroid hormones '
                                        'that are produced in the adrenal '
                                        'cortex.',
                          'parents': ['CHEBI:26764', 'CHEBI:50858']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 26,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}