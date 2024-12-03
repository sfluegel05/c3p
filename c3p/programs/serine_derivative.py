"""
Classifies: CHEBI:26649 serine derivative
"""
from rdkit import Chem

def is_serine_derivative(smiles: str):
    """
    Determines if a molecule is a serine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a serine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Serine core structure: N[C@H](CO)C(=O)O
    serine_smiles = "N[C@H](CO)C(=O)O"
    serine_mol = Chem.MolFromSmiles(serine_smiles)
    
    if not mol.HasSubstructMatch(serine_mol):
        return False, "Molecule does not contain the serine core structure"

    # Check for modifications at the amino group or carboxy group, or replacement of any hydrogen by a heteroatom
    # This is a simplified check, as a full check would require extensive substructure searching.
    
    amino_group = Chem.MolFromSmarts("[N]-[C@H](CO)C(=O)O")
    carboxy_group = Chem.MolFromSmarts("N[C@H](CO)C(=O)[O]")
    heteroatom_replacement = Chem.MolFromSmarts("N[C@H](CO)C(=O)O[*]")

    if mol.HasSubstructMatch(amino_group) or mol.HasSubstructMatch(carboxy_group) or mol.HasSubstructMatch(heteroatom_replacement):
        return True, "Molecule is a serine derivative"
    
    return False, "Molecule does not meet the criteria for a serine derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26649',
                          'name': 'serine derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of serine at the amino '
                                        'group or the carboxy group, or from '
                                        'the replacement of any hydrogen of '
                                        'serine by a heteroatom. The '
                                        'definition normally excludes peptides '
                                        'containing serine residues.',
                          'parents': ['CHEBI:83821']},
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
    'num_true_positives': 68,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}