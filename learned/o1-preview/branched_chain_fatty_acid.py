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
    from collections import deque

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid group (C(=O)O[H])
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    matches = mol.GetSubstructMatches(carboxylic_acid)
    if len(matches) == 0:
        return False, "No carboxylic acid group found"

    # Assume the first carboxyl carbon is the start
    carboxyl_idx = matches[0][0]

    # Remove carboxyl oxygens to focus on hydrocarbon chain
    mol_no_oxygen = Chem.RWMol(mol)
    for atom in mol_no_oxygen.GetAtoms():
        if atom.GetAtomicNum() == 8:
            atom.SetAtomicNum(0)  # Change oxygen to dummy atom
    mol_no_oxygen.UpdatePropertyCache()

    # Find the longest carbon chain starting from the carboxyl carbon
    def longest_chain_from_atom(mol, start_idx):
        visited = set()
        max_path = []

        def dfs(current_idx, path):
            nonlocal max_path
            visited.add(current_idx)
            path = path + [current_idx]
            if len(path) > len(max_path):
                max_path = path
            atom = mol.GetAtomWithIdx(current_idx)
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                nbr_atom = mol.GetAtomWithIdx(nbr_idx)
                if nbr_idx not in visited and nbr_atom.GetAtomicNum() == 6:
                    dfs(nbr_idx, path)
            visited.remove(current_idx)

        dfs(start_idx, [])
        return max_path

    main_chain = longest_chain_from_atom(mol_no_oxygen, carboxyl_idx)

    if len(main_chain) == 0:
        return False, "No hydrocarbon chain found"

    # Check for branching: atoms in the main chain with side chains
    branching_found = False
    for idx in main_chain:
        atom = mol_no_oxygen.GetAtomWithIdx(idx)
        neighbor_indices = [nbr.GetIdx() for nbr in atom.GetNeighbors()
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in main_chain]
        if neighbor_indices:
            branching_found = True
            break

    if branching_found:
        return True, "Branched-chain fatty acid detected"

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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}