"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: CHEBI:34084 aldopentose
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    An aldopentose is a pentose with a (potential) aldehyde group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldopentose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly 5 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 5:
        return False, f"Expected 5 carbons, found {c_count}"

    # Check for aldehyde group (either explicit or potential)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    # If no explicit aldehyde, check for potential aldehyde (e.g., in cyclic form)
    if not aldehyde_matches:
        # Look for a carbon with a single bond to oxygen (potential aldehyde in cyclic form)
        potential_aldehyde_pattern = Chem.MolFromSmarts("[CX4][OX2]")
        potential_aldehyde_matches = mol.GetSubstructMatches(potential_aldehyde_pattern)
        if not potential_aldehyde_matches:
            return False, "No aldehyde or potential aldehyde group found"

    # Check for multiple hydroxyl groups (at least 3)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if hydroxyl_count < 3:
        return False, f"Expected at least 3 hydroxyl groups, found {hydroxyl_count}"

    # Check molecular weight (should be around 150 g/mol for aldopentoses)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 200:
        return False, f"Molecular weight {mol_wt:.2f} is outside expected range for aldopentoses"

    return True, "Contains 5 carbons, a (potential) aldehyde group, and multiple hydroxyl groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:34084',
                          'name': 'aldopentose',
                          'definition': 'A pentose with a (potential) aldehyde '
                                        'group at one end.',
                          'parents': ['CHEBI:46983', 'CHEBI:26561']},
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