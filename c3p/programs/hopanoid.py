"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:36243 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the hopane skeleton pattern
    hopane_pattern = Chem.MolFromSmarts(
        "[C@H]1CC[C@@]2(C)[C@H]1CC[C@]1(C)[C@@H]2CC[C@@H]2[C@@]3(C)CCCC(C)(C)[C@@H]3CC[C@@]12C"
    )
    
    # Check if the molecule contains the hopane skeleton
    if not mol.HasSubstructMatch(hopane_pattern):
        return False, "No hopane skeleton found"

    # Check the number of rings to ensure it's a pentacyclic structure
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 5:
        return False, f"Found {n_rings} rings, need at least 5 for hopanoid"

    # Check molecular weight - hopanoids typically have a high molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for hopanoid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 30:
        return False, "Too few carbons for hopanoid"
    if o_count > 10:
        return False, "Too many oxygens for typical hopanoid"

    return True, "Contains hopane skeleton with appropriate ring structure and molecular weight"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36243',
                          'name': 'hopanoid',
                          'definition': 'A triterpenoid based on a hopane skeleton.',
                          'parents': ['CHEBI:25805', 'CHEBI:36243']},
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