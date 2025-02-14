"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:18148 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a 4-(2-aminoethyl)pyrocatechol (4-(2-aminoethyl)benzene-1,2-diol)
    or a derivative formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pyrocatechol (benzene-1,2-diol) core
    pyrocatechol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1")
    if not mol.HasSubstructMatch(pyrocatechol_pattern):
        return False, "No pyrocatechol core found"

    # Look for 2-aminoethyl side chain
    ethylamine_pattern = Chem.MolFromSmarts("CCN")
    if not mol.HasSubstructMatch(ethylamine_pattern):
        return False, "No 2-aminoethyl side chain found"

    # Check connectivity between pyrocatechol and ethylamine
    connected = AllChem.GetDistanceMatrix(mol).max() < 10  # Max topological distance
    if not connected:
        return False, "Pyrocatechol and ethylamine not connected"

    # Check for alternative substitutions
    if mol.GetSubstructMatches(pyrocatechol_pattern)[0][0].GetDegree() > 3:
        return True, "Substituted catecholamine derivative"

    return True, "Contains pyrocatechol core with 2-aminoethyl side chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:18148',
        'name': 'catecholamine',
        'definition': '4-(2-Aminoethyl)pyrocatechol [4-(2-aminoethyl)benzene-1,2-diol] and derivatives formed by substitution.',
        'parents': ['CHEBI:18669', 'CHEBI:24506', 'CHEBI:50113']
    },
    'config': {'llm_model_name': 'lbl/claude-sonnet', 'f1_threshold': 0.8, 'max_attempts': 5, 'max_positive_instances': None, 'max_positive_to_test': None, 'max_negative_to_test': None, 'max_positive_in_prompt': 50, 'max_negative_in_prompt': 20, 'max_instances_in_prompt': 100, 'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 144,
    'num_false_positives': 1,
    'num_true_negatives': 182418,
    'num_false_negatives': 21,
    'num_negatives': None,
    'precision': 0.9930555555555556,
    'recall': 0.8728323699421965,
    'f1': 0.9286389394717246,
    'accuracy': 0.9998880047278937
}