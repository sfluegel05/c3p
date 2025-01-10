"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: CHEBI:36044 lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for hydroperoxy group (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2][OX1]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"

    # Check for lipid characteristics (long carbon chain and carboxylic acid/ester)
    # Look for at least 10 carbons in a chain
    long_chain_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long carbon chain found"

    # Check for carboxylic acid or ester group
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid or ester group found"

    # Check molecular weight - lipids typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for lipid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for lipid"
    if o_count < 3:
        return False, "Too few oxygens for lipid hydroperoxide"

    return True, "Contains hydroperoxy group and lipid characteristics"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:36044',
        'name': 'lipid hydroperoxide',
        'definition': 'Any lipid carrying one or more hydroperoxy substituents.',
        'parents': ['CHEBI:36044', 'CHEBI:36044']
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