"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: CHEBI:57540 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate group ([O-]C=O)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Check for a long carbon chain (at least 4 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found"

    # Count the number of carbons in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
        return False, "Too few carbons to be a fatty acid anion"

    # Check molecular weight - fatty acid anions typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for fatty acid anion"

    return True, "Contains a carboxylate group and a long carbon chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:57540',
        'name': 'fatty acid anion',
        'definition': 'The conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group of the corresponding fatty acid.',
        'parents': ['CHEBI:57539', 'CHEBI:24834']
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