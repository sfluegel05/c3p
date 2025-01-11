"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: CHEBI:36235 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid core conjugated to a functional group that adds hydrophilicity or charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the bile acid core pattern (steroid-like structure with hydroxyl groups)
    bile_acid_core_pattern = Chem.MolFromSmarts("[C@H]1[C@@]2([C@H]([C@@H]3[C@]([C@@H]([C@H]4[C@](CC3)(C[C@H](O)CC4)[H])C)[H])CC[C@@]2([C@@H](C1)O)[H])")
    if not mol.HasSubstructMatch(bile_acid_core_pattern):
        return False, "No bile acid core structure found"

    # Define patterns for conjugated groups
    conjugated_groups = [
        Chem.MolFromSmarts("[NX3][CX3](=[OX1])"),  # Amino acid (e.g., glycine, taurine)
        Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])([OX2])"),  # Sulfate
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])"),  # Glucuronic acid
        Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])[CX6H1]([OX2H1])"),  # Glucose
        Chem.MolFromSmarts("[CX3](=[OX1])[NX3][CX3](=[OX1])"),  # Coenzyme A-like structure
    ]

    # Check for presence of conjugated groups
    has_conjugated_group = any(mol.HasSubstructMatch(pattern) for pattern in conjugated_groups)
    if not has_conjugated_group:
        return False, "No conjugated group found"

    # Check for hydrophilicity or charge
    # Calculate the number of polar atoms (O, N, S) and formal charges
    polar_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8, 16])
    formal_charges = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())

    if polar_atoms < 3 and formal_charges == 0:
        return False, "Insufficient hydrophilicity or charge"

    return True, "Contains bile acid core with conjugated hydrophilic or charged group"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36235',
        'name': 'bile acid conjugate',
        'definition': 'Any bile acid conjugated to a functional group that gives additional hydrophilicity or charge to the molecule.',
        'parents': ['CHEBI:36234', 'CHEBI:36236']
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
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}