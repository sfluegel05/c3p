"""
Classifies: CHEBI:26395 purine nucleotide
"""
from rdkit import Chem

def is_purine_nucleotide(smiles: str):
    """
    Determines if a molecule is a purine nucleotide (a nucleotide that has a purine nucleobase).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a purine nucleotide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    purine_smarts = Chem.MolFromSmarts('c1ncnc2n(cnc12)')
    if not mol.HasSubstructMatch(purine_smarts):
        return False, "No purine nucleobase found"

    phosphate_smarts = Chem.MolFromSmarts('P(=O)(O)O')
    if not mol.HasSubstructMatch(phosphate_smarts):
        return False, "No phosphate group found"

    ribose_smarts = Chem.MolFromSmarts('O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O')
    if not mol.HasSubstructMatch(ribose_smarts):
        return False, "No ribose sugar found"

    return True, "Molecule is a purine nucleotide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26395',
                          'name': 'purine nucleotide',
                          'definition': 'Any nucleotide that has a purine '
                                        'nucleobase.',
                          'parents': ['CHEBI:26401', 'CHEBI:36976']},
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
    'num_false_negatives': 24,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}