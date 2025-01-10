"""
Classifies: CHEBI:75769 B vitamin
"""
"""
Classifies: CHEBI:33284 B vitamin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_B_vitamin(smiles: str):
    """
    Determines if a molecule is a B vitamin based on its SMILES string.
    B vitamins are a group of water-soluble vitamins with specific structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a B vitamin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for nitrogen-containing heterocycles (common in B vitamins)
    heterocycle_patterns = [
        "[n]1ccccc1",  # Pyridine
        "[n]1ccncc1",  # Pyrimidine
        "[n]1cncc1",   # Imidazole
        "[n]1ccnc1",   # Pyrazole
        "[n]1cc[nH]c1" # Pyrrole
    ]
    has_heterocycle = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in heterocycle_patterns)
    if not has_heterocycle:
        return False, "No nitrogen-containing heterocycle found"

    # Check for functional groups common in B vitamins
    functional_group_patterns = [
        "[CX3](=O)[OX2H1]",  # Carboxylic acid
        "[NX3][CX4]",        # Amine
        "[PX4](=O)([OX2H1])",# Phosphate
        "[CX3](=O)[OX1H0-]", # Carboxylate
        "[CX3](=O)[NX3]",    # Amide
        "[SX2](=O)(=O)[OX2H1]" # Sulfonate
    ]
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in functional_group_patterns)
    if not has_functional_group:
        return False, "No functional group common in B vitamins found"

    # Check for water solubility (B vitamins are water-soluble)
    logP = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
    if logP > 1.5:  # Arbitrary threshold for water solubility
        return False, "Molecule is likely not water-soluble"

    # Check molecular weight (B vitamins typically have MW < 1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, "Molecular weight too high for a B vitamin"

    return True, "Contains nitrogen-containing heterocycle and functional groups common in B vitamins, and is water-soluble"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33284',
                          'name': 'B vitamin',
                          'definition': 'Any member of the group of eight water-soluble vitamins originally thought to be a single compound (vitamin B) that play important roles in cell metabolism. The group comprises of vitamin B1, B2, B3, B5, B6, B7, B9, and B12 (Around 20 other compounds were once thought to be B vitamins but are no longer classified as such).',
                          'parents': ['CHEBI:33280', 'CHEBI:33281']},
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