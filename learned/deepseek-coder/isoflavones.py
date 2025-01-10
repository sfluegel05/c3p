"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: CHEBI:50753 isoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core isoflavone pattern: 3-aryl-1-benzopyran-4-one
    # The pattern captures the essential features of the isoflavone structure:
    # - A benzopyran-4-one core (C1=CC(=O)C2=C(O1)C=CC=C2)
    # - An aryl group attached at the 3-position of the benzopyran-4-one core
    isoflavone_pattern = Chem.MolFromSmarts("[c:1]1[c:2][c:3][c:4][c:5][c:6]1-[c:7]2[c:8][c:9][c:10][c:11][c:12]2-[O:13]-[C:14]=[O:15]")
    
    # Check if the molecule matches the isoflavone pattern
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No 3-aryl-1-benzopyran-4-one skeleton found"

    # Verify the presence of an aryl group at the 3-position
    aryl_pattern = Chem.MolFromSmarts("[c:1]1[c:2][c:3][c:4][c:5][c:6]1")
    aryl_matches = mol.GetSubstructMatches(aryl_pattern)
    if len(aryl_matches) < 1:
        return False, "No aryl group found at the 3-position"

    # Check for the benzopyran-4-one core
    benzopyranone_pattern = Chem.MolFromSmarts("[c:1]1[c:2][c:3][c:4][c:5][c:6]1-[O:7]-[C:8]=[O:9]")
    if not mol.HasSubstructMatch(benzopyranone_pattern):
        return False, "No benzopyran-4-one core found"

    # Count the number of aromatic rings to ensure the presence of the aryl group
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if n_aromatic_rings < 2:
        return False, "Not enough aromatic rings for isoflavone structure"

    # Check molecular weight - isoflavones typically have MW > 200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for isoflavone"

    return True, "Contains 3-aryl-1-benzopyran-4-one skeleton"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50753',
                          'name': 'isoflavones',
                          'definition': 'Any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton and its substituted derivatives.',
                          'parents': ['CHEBI:50753', 'CHEBI:50753']},
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