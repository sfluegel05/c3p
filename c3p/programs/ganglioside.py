"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for key components of gangliosides
    ceramide_smarts = '[C@H](NC=O)[C@H](O)C'
    oligosaccharide_smarts = '[C@H]1O[C@H](CO)[C@H](O)[C@@H](O)[C@H]1O'
    sialic_acid_smarts = '[C@H](O)C(=O)[O-]'

    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)
    oligosaccharide_pattern = Chem.MolFromSmarts(oligosaccharide_smarts)
    sialic_acid_pattern = Chem.MolFromSmarts(sialic_acid_smarts)

    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide structure found"

    if not mol.HasSubstructMatch(oligosaccharide_pattern):
        return False, "No oligosaccharide structure found"

    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No sialic acid structure found"

    return True, "Ganglioside structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28892',
                          'name': 'ganglioside',
                          'definition': 'A molecule composed of a '
                                        'glycosphingolipid (ceramide and '
                                        'oligosaccharide) with one or more '
                                        'sialic acids linked on the sugar '
                                        'chain.',
                          'parents': [   'CHEBI:17761',
                                         'CHEBI:231691',
                                         'CHEBI:36526']},
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
    'num_false_negatives': 29,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}