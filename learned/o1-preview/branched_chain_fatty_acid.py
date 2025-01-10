"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    
    A branched-chain fatty acid is a fatty acid in which the hydrocarbon chain has one or more
    alkyl substituents (branches). The fatty acyl chain is usually saturated and the substituent
    is a methyl group; however, unsaturated BCFAs are found in marine animals, and branches
    other than methyl are found in microbial lipids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (C(=O)O[H])
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Get the carboxyl carbon atom
    carboxyl_carbon_idx = carboxy_matches[0][0]
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)

    # Find the alpha carbon (carbon adjacent to the carboxyl carbon)
    alpha_carbons = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if not alpha_carbons:
        return False, "No alpha carbon found"
    alpha_carbon = alpha_carbons[0]

    # Traverse the hydrocarbon chain starting from the alpha carbon
    visited = set()
    stack = [(alpha_carbon, carboxyl_carbon)]  # (current_atom, previous_atom)

    while stack:
        current_atom, prev_atom = stack.pop()
        visited.add(current_atom.GetIdx())
        
        # Get all carbon neighbors excluding the previous atom
        neighbors = [
            nbr for nbr in current_atom.GetNeighbors()
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()
        ]
        
        # Check if the current carbon atom is a branching point
        if len(neighbors) > 1:
            return True, "Branched chain found in hydrocarbon chain"
        
        # Add unvisited carbon neighbors to the stack
        for nbr in neighbors:
            if nbr.GetIdx() not in visited:
                stack.append((nbr, current_atom))

    # If no branching points are found, it's a linear fatty acid
    return False, "No branching in hydrocarbon chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25941',
        'name': 'branched-chain fatty acid',
        'definition': 'Any fatty acid in which the parent hydrocarbon chain has one or more alkyl substituents; a common component in animal and bacterial lipids. The fatty acyl chain is usually saturated and the substituent a methyl group; however, unsaturated BCFAs are found in marine animals, and branches other than methyl are found in microbial lipids.',
        'parents': ['CHEBI:15740', 'CHEBI:25359']
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