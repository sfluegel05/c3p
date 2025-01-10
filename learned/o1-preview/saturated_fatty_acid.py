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
    
    # Remove isotopes and explicit hydrogens for substructure matching
    mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H0]')
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
    
    # Count number of carbons in the main chain
    # Extract the carboxylic acid carbon
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"
    else:
        carboxylic_acid_carbon_idx = carboxylic_acid_matches[0][0]

    # Use a BFS to find the longest carbon chain starting from the carboxylic acid carbon
    visited = set()
    max_chain_length = 0
    from collections import deque
    queue = deque([(carboxylic_acid_carbon_idx, 0)])
    while queue:
        atom_idx, length = queue.popleft()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            max_chain_length = max(max_chain_length, length)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    queue.append((neighbor_idx, length + 1))
    
    # A fatty acid typically has at least 3 carbons in the chain excluding the carboxyl carbon
    if max_chain_length < 3:
        return False, f"Carbon chain too short for fatty acid (found {max_chain_length} carbons)"
    
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}