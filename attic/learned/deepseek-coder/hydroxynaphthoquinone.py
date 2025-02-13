"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:XXXXX hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone with at least one hydroxy group attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general naphthoquinone core pattern
    naphthoquinone_pattern = Chem.MolFromSmarts("c1ccc2C(=O)C=CC(=O)c2c1 |c1:1,2,3,4,5,6,7,8|")
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        # Try another possible naphthoquinone pattern (1,2-naphthoquinone)
        naphthoquinone_pattern = Chem.MolFromSmarts("c1ccc2C(=O)C(=O)C=Cc2c1 |c1:1,2,3,4,5,6,7,8|")
        if not mol.HasSubstructMatch(naphthoquinone_pattern):
            return False, "No naphthoquinone core found"

    # Define the hydroxy group pattern
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) == 0:
        return False, "No hydroxy group found"

    # Check if at least one hydroxy group is directly attached to the naphthoquinone core
    naphthoquinone_atoms = mol.GetSubstructMatch(naphthoquinone_pattern)
    for match in hydroxy_matches:
        hydroxy_atom = match[0]
        for neighbor in mol.GetAtomWithIdx(hydroxy_atom).GetNeighbors():
            if neighbor.GetIdx() in naphthoquinone_atoms:
                return True, "Contains a naphthoquinone core with at least one hydroxy group attached"

    return False, "Hydroxy group not directly attached to the naphthoquinone core"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXXX',
                          'name': 'hydroxynaphthoquinone',
                          'definition': 'Any naphthoquinone in which the naphthoquinone moiety is substituted by at least one hydroxy group.',
                          'parents': ['CHEBI:XXXXX', 'CHEBI:XXXXX']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}