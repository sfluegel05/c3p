"""
Classifies: CHEBI:47003 leucine derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_leucine_derivative(smiles: str):
    """
    Determines if a molecule is a leucine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leucine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the leucine structure
    leucine_smiles = "CC(C)CC(C(=O)O)N"
    leucine = Chem.MolFromSmiles(leucine_smiles)
    if leucine is None:
        return False, "Failed to create leucine molecule"

    # Perform substructure search to find leucine core
    if not mol.HasSubstructMatch(leucine):
        return False, "Leucine core structure not found"

    # Check for modifications at the amino group or carboxy group
    amino_group = Chem.MolFromSmarts("[NX3H2,NX3H1,NX3H0]")
    carboxy_group = Chem.MolFromSmarts("C(=O)[O-,OH]")
    modified_amino = mol.HasSubstructMatch(amino_group)
    modified_carboxy = mol.HasSubstructMatch(carboxy_group)

    if modified_amino or modified_carboxy:
        return True, "Leucine derivative with modified amino or carboxy group"

    # Check for replacement of any hydrogen by a heteroatom
    heteroatoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6]]
    if heteroatoms:
        return True, "Leucine derivative with heteroatom substitution"

    return False, "No modifications found that classify as leucine derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47003',
                          'name': 'leucine derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of leucine at the amino '
                                        'group or the carboxy group, or from '
                                        'the replacement of any hydrogen of '
                                        'leucine by a heteroatom. The '
                                        'definition normally excludes peptides '
                                        'containing leucine residues.',
                          'parents': ['CHEBI:83821']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 8,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 4,
    'precision': 1.0,
    'recall': 0.6666666666666666,
    'f1': 0.8,
    'accuracy': None}