"""
Classifies: CHEBI:35436 D-glucoside
"""
"""
Classifies: D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is any glucoside in which the glycosyl group is derived from D-glucose.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the D-glucoside pattern with correct stereochemistry
    # D-glucose unit connected via O-glycosidic bond at C1
    d_glucoside_smarts = """
    [C@H]1([O][#6])
    [O][C@@H]([C@H]([C@@H]([C@H]1O)O)O)O
    """
    pattern = Chem.MolFromSmarts(d_glucoside_smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern for D-glucoside"
    
    # Check for D-glucoside substructure with correct stereochemistry
    matches = mol.GetSubstructMatches(pattern, useChirality=True)
    if not matches:
        return False, "No D-glucoside moiety with correct stereochemistry found"
    
    return True, "Contains D-glucoside moiety with correct stereochemistry"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:24269',
        'name': 'D-glucoside',
        'definition': 'Any glucoside in which the glycosyl group is derived from D-glucose.',
        'parents': ['CHEBI:48226', 'CHEBI:26648']
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Classification metrics would be filled in after testing the function
}