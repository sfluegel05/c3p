"""
Classifies: CHEBI:62761 tyrosine derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_tyrosine_derivative(smiles: str):
    """
    Determines if a molecule is a tyrosine derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tyrosine derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core structure of tyrosine
    tyrosine_smiles = "O=C(O)[C@@H](N)Cc1ccc(O)cc1"
    tyrosine_mol = Chem.MolFromSmiles(tyrosine_smiles)
    if tyrosine_mol is None:
        return False, "Invalid core tyrosine structure"

    # Check if the molecule contains the tyrosine core structure
    if not mol.HasSubstructMatch(tyrosine_mol):
        return False, "Does not contain core tyrosine structure"

    # Check for modifications at the amino group, carboxy group, or phenyl ring
    amino_group = Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")
    carboxy_group = Chem.MolFromSmarts("C(=O)[O-,OH]")
    phenyl_ring = Chem.MolFromSmarts("c1ccccc1")
    
    if mol.HasSubstructMatch(amino_group) or mol.HasSubstructMatch(carboxy_group) or mol.HasSubstructMatch(phenyl_ring):
        return True, "Contains modifications at the amino group, carboxy group, or phenyl ring"

    # Check for replacement of any hydrogen by a heteroatom
    heteroatom_substitution = Chem.MolFromSmarts("[!#1]")
    if mol.HasSubstructMatch(heteroatom_substitution):
        return True, "Contains heteroatom substitution"

    return False, "No relevant modifications found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62761',
                          'name': 'tyrosine derivative',
                          'definition': 'An amino acid derivative resulting '
                                        'from reaction of tyrosine at the '
                                        'amino group or the carboxy group, any '
                                        'substitution of phenyl hydrogens, or '
                                        'from the replacement of any hydrogen '
                                        'of tyrosine by a heteroatom. The '
                                        'definition normally excludes peptides '
                                        'containing tyrosine residues.',
                          'parents': ['CHEBI:25985', 'CHEBI:83821']},
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
    'num_true_positives': 13,
    'num_false_positives': 3,
    'num_true_negatives': 12,
    'num_false_negatives': 2,
    'precision': 0.8125,
    'recall': 0.8666666666666667,
    'f1': 0.8387096774193549,
    'accuracy': None}