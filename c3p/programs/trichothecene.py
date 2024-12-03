"""
Classifies: CHEBI:55517 trichothecene
"""
from rdkit import Chem

def is_trichothecene(smiles: str):
    """
    Determines if a molecule is a trichothecene.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a trichothecene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a 12,13-epoxy ring
    epoxy_pattern = Chem.MolFromSmarts('C1C(O1)C')
    if not mol.HasSubstructMatch(epoxy_pattern):
        return False, "No 12,13-epoxy ring found"

    # Check for the presence of hydroxy or acetyl groups
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    acetyl_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(hydroxy_pattern) and not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No hydroxy or acetyl groups found"

    # Check for sesquiterpene skeleton
    sesquiterpene_pattern = Chem.MolFromSmarts('C1CC2CCCC3CCCC(C1)C23')
    if not mol.HasSubstructMatch(sesquiterpene_pattern):
        return False, "No sesquiterpene skeleton found"

    return True, "Molecule is a trichothecene"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:55517',
                          'name': 'trichothecene',
                          'definition': 'Any one of a large family of '
                                        'chemically related mycotoxins with a '
                                        'structure based on a sesquiterpene '
                                        'skeleton. The most important '
                                        'structural features causing the '
                                        'biological activities of '
                                        'trichothecenes are a 12,13-epoxy '
                                        'ring, the presence of hydroxy or '
                                        'acetyl groups at appropriate '
                                        'positions on the trichothecene '
                                        'nucleus and the structure and '
                                        'position of the side-chain.',
                          'parents': [   'CHEBI:26658',
                                         'CHEBI:32955',
                                         'CHEBI:38166']},
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
    'num_true_negatives': 16,
    'num_false_negatives': 16,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}