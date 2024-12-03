"""
Classifies: CHEBI:51270 tetracenes
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_tetracenes(smiles: str):
    """
    Determines if a molecule is a tetracene (compounds containing a tetracene skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetracene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Tetracene consists of four linearly-fused benzene rings
    tetracene_smarts = 'c1ccc2c(c1)ccc3c2ccc4c(c3)cccc4'
    tetracene_pattern = Chem.MolFromSmarts(tetracene_smarts)

    if mol.HasSubstructMatch(tetracene_pattern):
        return True, "Molecule contains a tetracene skeleton"
    else:
        return False, "Molecule does not contain a tetracene skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51270',
                          'name': 'tetracenes',
                          'definition': 'Compounds containing a tetracene '
                                        'skeleton.',
                          'parents': ['CHEBI:51269']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 20,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}