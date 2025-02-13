"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: CHEBI:15841 saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Check for carboxylic acid group (-C(=O)OH or deprotonated -C(=O)O-)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O,H,-]')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for carbon-carbon double bonds
    if mol.HasSubstructMatch(Chem.MolFromSmarts('C=C')):
        return False, "Contains carbon-carbon double bonds (unsaturation)"

    # Check for carbon-carbon triple bonds
    if mol.HasSubstructMatch(Chem.MolFromSmarts('C#C')):
        return False, "Contains carbon-carbon triple bonds (unsaturation)"
        
    # Check for ring structures
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains ring structures, not a fatty acid"
    
    # Find the carboxylic acid carbon atom
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxylic_acid_carbons = [match[0] for match in carboxylic_acid_matches]
    
    if not carboxylic_acid_carbons:
        return False, "No carboxylic acid group found"
    
    # Find the longest carbon chain attached to the carboxylic acid group
    max_chain_length = 0
    for carboxylic_carbon_idx in carboxylic_acid_carbons:
        chain_length = get_longest_aliphatic_chain(mol, carboxylic_carbon_idx)
        max_chain_length = max(max_chain_length, chain_length)
    
    # A fatty acid should have at least 4 carbons in the chain excluding the carboxyl carbon
    if max_chain_length < 4:
        return False, f"Carbon chain too short for fatty acid (found {max_chain_length} carbons)"
    
    return True, "Molecule is a saturated fatty acid"

def get_longest_aliphatic_chain(mol, start_atom_idx):
    """
    Finds the length of the longest aliphatic carbon chain starting from the given atom.

    Args:
        mol: RDKit molecule object
        start_atom_idx: Index of the starting atom (carboxylic acid carbon)

    Returns:
        int: Length of the longest carbon chain
    """
    visited = set()
    max_length = 0
    stack = [(start_atom_idx, 0)]
    
    while stack:
        atom_idx, length = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        
        if atom.GetAtomicNum() != 6:
            continue  # Only consider carbon atoms
        
        max_length = max(max_length, length)
        
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            # Skip if bond is not single
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            if neighbor_idx not in visited:
                stack.append((neighbor_idx, length + 1))
                
    return max_length

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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}