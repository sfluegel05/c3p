"""
Classifies: CHEBI:23509 cysteine derivative
"""
from rdkit import Chem

def is_cysteine_derivative(smiles: str):
    """
    Determines if a molecule is a cysteine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cysteine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the cysteine core structure
    cysteine_core = Chem.MolFromSmiles("N[C@@H](CS)C(=O)O")
    if cysteine_core is None:
        return False, "Error in defining cysteine core structure"

    # Check if the molecule contains the cysteine core
    if not mol.HasSubstructMatch(cysteine_core):
        return False, "Molecule does not contain cysteine core structure"

    # Check for modifications at the amino group, carboxy group, or thiol group
    amino_modifications = Chem.MolFromSmarts("N[C@@H](CS)C(=O)O")
    carboxy_modifications = Chem.MolFromSmarts("N[C@@H](CS)C(=O)O")
    thiol_modifications = Chem.MolFromSmarts("N[C@@H](CS)C(=O)O")

    if mol.HasSubstructMatch(amino_modifications) or mol.HasSubstructMatch(carboxy_modifications) or mol.HasSubstructMatch(thiol_modifications):
        return True, "Molecule is a cysteine derivative with modifications at the amino, carboxy, or thiol group"

    # Check for replacement of any hydrogen by a heteroatom
    heteroatoms = ['N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I']
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in heteroatoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'H':
                    return True, "Molecule is a cysteine derivative with heteroatom substitution"

    return False, "Molecule does not meet the criteria for a cysteine derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23509',
                          'name': 'cysteine derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of cysteine at the '
                                        'amino group, carboxy group, or thiol '
                                        'group, or from the replacement of any '
                                        'hydrogen of cysteine by a heteroatom. '
                                        'The definition normally excludes '
                                        'peptides containing cysteine '
                                        'residues.',
                          'parents': ['CHEBI:83821']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 13,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.7692307692307693,
    'f1': 0.8695652173913044,
    'accuracy': None}