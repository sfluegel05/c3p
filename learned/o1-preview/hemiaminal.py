"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:51964 hemiaminal
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic compound where a carbon atom has both an amino group
    and a hydroxy group attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define hemiaminal SMARTS pattern
    # Carbon connected to OH group and nitrogen (single bonds)
    hemiaminal_pattern = Chem.MolFromSmarts("[C;!$(C=O);!$(C=N);!$(C#N)]([O;H1])([N;!$(N=*);!$(N#*)])")

    if hemiaminal_pattern is None:
        return False, "Invalid SMARTS pattern for hemiaminal"

    # Search for hemiaminal pattern in molecule
    if mol.HasSubstructMatch(hemiaminal_pattern):
        return True, "Contains hemiaminal functional group (carbon with attached OH and NH2/NHR/NR2)"
    else:
        return False, "No hemiaminal functional group found"

__metadata__ = {  
    'chemical_class': {   
        'id': 'CHEBI:51964',
        'name': 'hemiaminal',
        'definition': 'Any organic amino compound that has an amino group and a hydroxy group attached to the same carbon atom. Hemiaminals are intermediates in the formation of imines by addition of an amine to an aldehyde or ketone; those derived from primary amines are particularly unstable.',
        'parents': ['CHEBI:33709', 'CHEBI:49114']
    },
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
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}