"""
Classifies: CHEBI:50941 azaphilone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_azaphilone(smiles: str):
    """
    Determines if a molecule is an azaphilone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azaphilone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure patterns for azaphilone core structures
    isochromene_dione_pattern = Chem.MolFromSmarts('O=C1C=CC2=C(O1)C=CC(=O)C2')
    isoquinoline_dione_pattern = Chem.MolFromSmarts('O=C1C=CC2=C(O1)C=CC(=O)N2')

    if mol.HasSubstructMatch(isochromene_dione_pattern):
        return True, "Contains 6H-isochromene-6,8(7H)-dione skeleton"
    elif mol.HasSubstructMatch(isoquinoline_dione_pattern):
        return True, "Contains isoquinoline-6,8(2H,7H)-dione skeleton"
    else:
        return False, "Does not contain azaphilone core structure"

# Example usage:
# smiles = "O=C1C2=C(OC(=C1)C)C(=C(O)C(=C2O)CC3=C4C(=COC(=C4)C)[C@H](O)[C@](C3=O)(O)C)C"
# result, reason = is_azaphilone(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50941',
                          'name': 'azaphilone',
                          'definition': 'Any member of a family of natural '
                                        'products which contains a '
                                        '6H-isochromene-6,8(7H)-dione or an '
                                        'isoquinoline-6,8(2H,7H)-dione '
                                        'skeleton and its substituted '
                                        'derivatives thereof.',
                          'parents': ['CHEBI:24532']},
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
    'num_false_negatives': 20,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}