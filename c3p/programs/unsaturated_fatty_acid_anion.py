"""
Classifies: CHEBI:2580 unsaturated fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate anion group [C(=O)[O-]]
    carboxylate_anion = Chem.MolFromSmarts('[C](=O)[O-]')
    if not mol.HasSubstructMatch(carboxylate_anion):
        return False, "No carboxylate anion group found"

    # Check for at least one C-C double bond
    double_bond = Chem.MolFromSmarts('C=C')
    if not mol.HasSubstructMatch(double_bond):
        return False, "No C-C double bond found"

    # Check if it is a fatty acid (long carbon chain)
    carbon_chain = Chem.MolFromSmarts('CCCCCCCC')
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chain found"

    return True, "Molecule is an unsaturated fatty acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:2580',
                          'name': 'unsaturated fatty acid anion',
                          'definition': 'Any fatty acid anion containing at '
                                        'least one C-C unsaturated bond; '
                                        'formed by deprotonation of the '
                                        'carboxylic acid moiety.',
                          'parents': ['CHEBI:28868']},
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
    'num_true_positives': 30,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 31,
    'precision': 1.0,
    'recall': 0.4918032786885246,
    'f1': 0.6593406593406593,
    'accuracy': None}