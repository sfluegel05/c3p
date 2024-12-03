"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline or its hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline or its hydrogenated derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for beta-carboline core
    beta_carboline_smarts = "c1c2c3[nH]c4cccc(c4c1CCN2)C3"

    # Check if the molecule matches the beta-carboline core
    beta_carboline_core = Chem.MolFromSmarts(beta_carboline_smarts)
    if mol.HasSubstructMatch(beta_carboline_core):
        return True, "Contains beta-carboline core"

    # Check for hydrogenated derivatives (tetrahydro-beta-carboline)
    tetrahydro_beta_carboline_smarts = "C1CN2CCc3c[nH]c4cccc(c34)C2=C1"
    tetrahydro_beta_carboline_core = Chem.MolFromSmarts(tetrahydro_beta_carboline_smarts)
    if mol.HasSubstructMatch(tetrahydro_beta_carboline_core):
        return True, "Contains tetrahydro-beta-carboline core"

    return False, "Does not contain beta-carboline or its hydrogenated derivatives"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:60834',
                          'name': 'beta-carbolines',
                          'definition': 'Any pyridoindole containing a '
                                        'beta-carboline skeleton and their '
                                        'hydrogenated derivatives',
                          'parents': ['CHEBI:48888']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[01:42:52] SMILES Parse Error: syntax error while parsing: '
             'O=C\x01NC(CC(O)(C(=O)O)C(C)C)C(/C1=C(\\O)/C=C/C=C/CCC)=O\n'
             '[01:42:52] SMILES Parse Error: Failed parsing SMILES '
             "'O=C\x01NC(CC(O)(C(=O)O)C(C)C)C(/C1=C(\\O)/C=C/C=C/CCC)=O' for "
             'input: '
             "'O=C\x01NC(CC(O)(C(=O)O)C(C)C)C(/C1=C(\\O)/C=C/C=C/CCC)=O'\n",
    'stdout': '',
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 98,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}