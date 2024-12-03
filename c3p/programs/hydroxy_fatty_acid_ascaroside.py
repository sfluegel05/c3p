"""
Classifies: CHEBI:78799 hydroxy fatty acid ascaroside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid_ascaroside(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid ascaroside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid ascaroside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the ascaroside core structure
    ascaroside_core = Chem.MolFromSmarts('C[C@H](O[C@@H]1O[C@@H](C)[C@H](O)C[C@H]1O)')
    if not mol.HasSubstructMatch(ascaroside_core):
        return False, "Ascaroside core structure not found"

    # Check for the hydroxy fatty acid part
    hydroxy_fatty_acid = Chem.MolFromSmarts('CCCCCCCCCCCCCC[C@H](O)CC(O)=O')
    if mol.HasSubstructMatch(hydroxy_fatty_acid):
        return True, "Molecule is a hydroxy fatty acid ascaroside"

    # Check for alternative hydroxy fatty acid structures
    alternative_hydroxy_fatty_acid_patterns = [
        'CCCCCCCCCCCCCCC\\C=C\\C(O)=O',
        'CCCCCCCCCCCC(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O',
        'CCCCCC(O)=O',
        'CCCCCCCC\\C=C\\C(O)=O',
        'CCCCCCCCC(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O',
        'CCCC(O)=O',
        'OCCCCCCC(O)=O',
        'OCCCCCCC\\C=C\\C(O)=O',
        'OCCCCCCCCCCCCCC(O)=O',
        'CCCCCCCCCCCCCCCCCCCC/C=C(/C(=O)O)\\C',
        'CCCCC\\C=C\\C(O)=O',
        'OCCCCCCCC(O)=O',
        'CCCCCCCCCCCCCCC[C@@H](O)CC(O)=O'
    ]
    
    for pattern in alternative_hydroxy_fatty_acid_patterns:
        hydroxy_fatty_acid = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(hydroxy_fatty_acid):
            return True, "Molecule is a hydroxy fatty acid ascaroside"

    return False, "Hydroxy fatty acid structure not found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:78799',
                          'name': 'hydroxy fatty acid ascaroside',
                          'definition': 'A glycoside resulting from the '
                                        'reaction of the alcoholic hydroxy '
                                        'group of a hydroxy fatty acid with '
                                        'alpha-ascarylopyranose or its '
                                        'derivatives. They function as '
                                        'semiochemicals in the nematode '
                                        'Caenorhabditis elegans.',
                          'parents': ['CHEBI:61697', 'CHEBI:79202']},
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
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 5,
    'precision': 1.0,
    'recall': 0.6666666666666666,
    'f1': 0.8,
    'accuracy': None}