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
    This includes compounds like flavonoids, coumarins, lignins, and others with a C6-C3 backbone.

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

    # Define generalized SMARTS patterns for phenylpropanoid structures
    patterns = {
        'phenylpropane_core': Chem.MolFromSmarts('c1ccccc1-C-C-C'),  # Phenyl ring linked to 3-carbon chain
        'phenylpropene_core': Chem.MolFromSmarts('c1ccccc1-C=C-C'),  # Phenyl ring linked to propenyl chain
        'cinnamic_acid': Chem.MolFromSmarts('c1ccccc1-C=C-C(=O)[O,N]'),  # Cinnamic acid scaffold
        'coumarin_core': Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1'),  # Coumarin core
        'flavonoid_core': Chem.MolFromSmarts('c1cc(-c2coc3c(=O)cc(-c4ccccc4)cc3c2=O)ccc1'),  # Flavonoid skeleton
        'chalcone_core': Chem.MolFromSmarts('c1ccccc1-C=C-C(=O)-c2ccccc2'),  # Chalcone scaffold
        'lignin_precursor': Chem.MolFromSmarts('c1ccccc1-C-C-O'),  # Phenylpropanoid alcohols
        'stilbene_core': Chem.MolFromSmarts('c1ccccc1-C=C-c2ccccc2'),  # Stilbene scaffold
    }

    # Check for matches to any of the patterns
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name.replace('_', ' ')} substructure"

    # Check for phenylpropanoid ethers and esters
    phenylpropanoid_ether = Chem.MolFromSmarts('c1ccccc1-C-C-O')
    phenylpropanoid_ester = Chem.MolFromSmarts('c1ccccc1-C-C-C(=O)O')
    if mol.HasSubstructMatch(phenylpropanoid_ether):
        return True, "Contains phenylpropanoid ether substructure"
    if mol.HasSubstructMatch(phenylpropanoid_ester):
        return True, "Contains phenylpropanoid ester substructure"

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