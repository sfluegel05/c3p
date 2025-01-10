"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
"""
Classifies: CHEBI:32395 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid based on its SMILES string.
    A monounsaturated fatty acid has one double or triple bond in the carbon chain and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Count double and triple bonds in the carbon chain
    double_bonds = 0
    triple_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bonds += 1
        elif bond.GetBondType() == Chem.BondType.TRIPLE:
            triple_bonds += 1

    total_unsaturations = double_bonds + triple_bonds

    # Check for exactly one double or triple bond
    if total_unsaturations != 1:
        return False, f"Found {total_unsaturations} unsaturations, need exactly 1"

    # Check that the rest of the chain is singly bonded
    # We can verify this by ensuring that the number of rotatable bonds is consistent with a long carbon chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:  # Minimum number of rotatable bonds for a fatty acid
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - fatty acids typically >100 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for a fatty acid"

    return True, "Contains a carboxylic acid group and exactly one double or triple bond in the carbon chain"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32395',
                          'name': 'monounsaturated fatty acid',
                          'definition': 'Any fatty acid with one double or triple bond in the fatty acid chain and singly bonded carbon atoms in the rest of the chain. MUFAs have positive effects on the cardiovascular system, and in diabetes treatment.'},
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