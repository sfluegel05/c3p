"""
Classifies: CHEBI:59874 N-acyl-L-alpha-amino acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_N_acyl_L_alpha_amino_acid_anion(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate group [C(=O)[O-]]
    carboxylate = Chem.MolFromSmarts('C(=O)[O-]')
    if not mol.HasSubstructMatch(carboxylate):
        return False, "No carboxylate group found"

    # Check for alpha amino acid structure [C@H](N)
    alpha_amino_acid = Chem.MolFromSmarts('[C@H](N)')
    if not mol.HasSubstructMatch(alpha_amino_acid):
        return False, "No alpha amino acid structure found"

    # Check for N-acyl structure [N-C(=O)]
    n_acyl = Chem.MolFromSmarts('N-C(=O)')
    if not mol.HasSubstructMatch(n_acyl):
        return False, "No N-acyl structure found"

    return True, "Molecule is an N-acyl-L-alpha-amino acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59874',
                          'name': 'N-acyl-L-alpha-amino acid anion',
                          'definition': 'A carboxylic acid anion that is the '
                                        'conjugate base of an '
                                        'N-acyl-L-alpha-amino acid arising '
                                        'from deprotonation of the C-1 carboxy '
                                        'group.',
                          'parents': ['CHEBI:29067']},
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
    'num_true_positives': 11,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9166666666666666,
    'f1': 0.9565217391304348,
    'accuracy': None}