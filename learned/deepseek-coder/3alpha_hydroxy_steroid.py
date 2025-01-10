"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:17490 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid has a steroid backbone with a hydroxyl group at the 3-position in the alpha orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general steroid backbone pattern (four fused rings)
    steroid_backbone_pattern = Chem.MolFromSmarts("[C@H]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2)[C@@H]4CCCC[C@]34C")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        # Try a more general pattern without specific stereochemistry
        steroid_backbone_pattern = Chem.MolFromSmarts("C12CCC3C(C1CCC2)CCC4=CCCCC34")
        if not mol.HasSubstructMatch(steroid_backbone_pattern):
            return False, "No steroid backbone found"

    # Define the 3alpha-hydroxy pattern (hydroxyl group at position 3 in alpha orientation)
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@H]1([C@@H](O)[C@H]2CC[C@H]3[C@H]([C@@H]1CC[C@@H]2)[C@@H]4CCCC[C@]34C)")
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        # Try a more general pattern without specific stereochemistry
        alpha_hydroxy_pattern = Chem.MolFromSmarts("C1(C(O)C2CCC3C(C1CCC2)CCC4=CCCCC34)")
        if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
            return False, "No 3alpha-hydroxy group found"

    return True, "Contains a steroid backbone with a 3alpha-hydroxy group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17490',
                          'name': '3alpha-hydroxy steroid',
                          'definition': 'A 3-hydroxy steroid in which the 3-hydroxy substituent is in the alpha-position.',
                          'parents': ['CHEBI:47857', 'CHEBI:47856']},
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