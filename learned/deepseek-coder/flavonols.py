"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:00000 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone with a hydroxyl group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavonol pattern: a flavone backbone with a hydroxyl group at position 3
    flavonol_pattern = Chem.MolFromSmarts("[O;H1]c1c(=O)c2c(O)cc(O)cc2oc1")
    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "No flavonol backbone with hydroxyl at position 3 found"

    # Check for the presence of the chromen-4-one ring (flavone backbone)
    chromen_pattern = Chem.MolFromSmarts("O=C1C=CC(=O)C2=C1C=CC=C2")
    if not mol.HasSubstructMatch(chromen_pattern):
        return False, "No chromen-4-one ring (flavone backbone) found"

    # Check for at least one hydroxyl group on the benzene rings
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Not enough hydroxyl groups on the benzene rings"

    # Check for the presence of the hydroxyl group at position 3
    position_3_hydroxyl_pattern = Chem.MolFromSmarts("[O;H1]c1c(=O)c2c(O)cc(O)cc2oc1")
    if not mol.HasSubstructMatch(position_3_hydroxyl_pattern):
        return False, "No hydroxyl group at position 3 of the heterocyclic ring"

    return True, "Contains a flavone backbone with a hydroxyl group at position 3 and additional hydroxyl groups on the benzene rings"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:00000',
                          'name': 'flavonols',
                          'definition': 'Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic ring is replaced by a hydroxy group.',
                          'parents': ['CHEBI:00000', 'CHEBI:00000']},
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