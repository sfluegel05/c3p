"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: CHEBI:28801 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid has a steroid core with a hydroxyl group at the 3rd position in the beta orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general steroid core pattern (four fused rings)
    steroid_core_pattern = Chem.MolFromSmarts("[C@H]12CC[C@H]3[C@@H](CC[C@H]4[C@@H](CC[C@H]34)C)[C@@]1(CC2)C")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        # Try a more general pattern without specific stereochemistry
        steroid_core_pattern = Chem.MolFromSmarts("C12CC[C@H]3C[C@@H](CC[C@H]4C[C@@H](CC[C@H]34)C)[C@@]1(CC2)C")
        if not mol.HasSubstructMatch(steroid_core_pattern):
            return False, "No steroid core structure found"

    # Define the 3beta-hydroxy pattern (hydroxyl at position 3 in beta orientation)
    beta_hydroxy_pattern = Chem.MolFromSmarts("[C@H]1([C@@H](O)CC[C@H]2[C@@H]1CC[C@H]3[C@@H]2CC[C@H]4[C@@H]3CC[C@H]4C)C")
    if not mol.HasSubstructMatch(beta_hydroxy_pattern):
        # Try a more general pattern without specific stereochemistry
        beta_hydroxy_pattern = Chem.MolFromSmarts("C1(C(O)CCC2C1CCC3C2CCC4C3CCC4C)C")
        if not mol.HasSubstructMatch(beta_hydroxy_pattern):
            return False, "No 3beta-hydroxy group found"

    # Check for the presence of a hydroxyl group at the 3rd position
    hydroxyl_at_3 = Chem.MolFromSmarts("[C@H]1([C@@H](O)CC[C@H]2[C@@H]1CC[C@H]3[C@@H]2CC[C@H]4[C@@H]3CC[C@H]4C)C")
    if not mol.HasSubstructMatch(hydroxyl_at_3):
        # Try a more general pattern without specific stereochemistry
        hydroxyl_at_3 = Chem.MolFromSmarts("C1(C(O)CCC2C1CCC3C2CCC4C3CCC4C)C")
        if not mol.HasSubstructMatch(hydroxyl_at_3):
            return False, "No hydroxyl group at the 3rd position"

    return True, "Contains a steroid core with a 3beta-hydroxy group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28801',
                          'name': '3beta-hydroxy steroid',
                          'definition': 'A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.',
                          'parents': ['CHEBI:28800', 'CHEBI:28802']},
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