"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:15841 saturated fatty acid
"""
from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    A saturated fatty acid is a fatty acid containing no carbon-to-carbon multiple bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H1]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for carbon-carbon multiple bonds (unsaturation)
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Check if both atoms are carbon
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
            # Check for double or triple bonds
            if bond.GetBondTypeAsDouble() > 1.0:
                return False, "Contains carbon-carbon multiple bonds (unsaturation)"
    
    # Optional: Check for sufficient carbon chain length (e.g., at least 3 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return False, f"Carbon chain too short for fatty acid (found {c_count} carbons)"
    
    return True, "Molecule is a saturated fatty acid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15841',
        'name': 'saturated fatty acid',
        'definition': 'Any fatty acid containing no carbon to carbon multiple bonds.',
        'parents': ['CHEBI:35366', 'CHEBI:26559']
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