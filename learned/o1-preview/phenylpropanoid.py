"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is any organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for phenylpropanoid structures
    patterns = {
        'phenylpropane': Chem.MolFromSmarts('c1ccccc1CCC'),  # Phenylpropane skeleton
        'cinnamic_acid': Chem.MolFromSmarts('c1ccccc1C=CC(=O)O'),  # Cinnamic acid scaffold
        'coumarin': Chem.MolFromSmarts('O=c1oc2ccccc2c(=O)c1'),  # Coumarin core
        'flavonoid': Chem.MolFromSmarts('c1cc(-c2ccc(O)cc2)c(=O)c3coc(=O)cc13'),  # Flavonoid skeleton
        'lignin_precursor': Chem.MolFromSmarts('c1ccccc1CCO'),  # Benzene ring with propanol side chain
    }

    # Check for matches to any of the patterns
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} substructure"

    return False, "Does not contain phenylpropanoid substructure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'phenylpropanoid',
        'definition': 'Any organic aromatic compound with a structure based on a phenylpropane skeleton. The class includes naturally occurring phenylpropanoid esters, flavonoids, anthocyanins, coumarins and many small phenolic molecules as well as their semi-synthetic and synthetic analogues. Phenylpropanoids are also precursors of lignin.',
        'parents': []},
    'config': {   
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}