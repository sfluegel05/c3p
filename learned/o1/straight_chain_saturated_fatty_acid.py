"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
"""
Classifies: CHEBI:26666 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    A straight-chain saturated fatty acid is an unbranched, fully saturated carbon chain ending with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Identify the carboxylic acid carbon
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    carboxylic_carbons = [match[0] for match in matches]  # First carbon in the pattern
    
    # Check that molecule is acyclic (no rings)
    if Chem.GetSSSR(mol) > 0:
        return False, "Molecule contains rings"
    
    # Check for unsaturation (double/triple bonds excluding carboxylic acid group)
    unsaturated_bond = False
    for bond in mol.GetBonds():
        bond_type = bond.GetBondType()
        if bond_type in (Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE):
            begin_atom = bond.GetBeginAtomIdx()
            end_atom = bond.GetEndAtomIdx()
            # Allow double bond in the carboxylic acid group (C=O)
            if bond_type == Chem.rdchem.BondType.DOUBLE and (
                (begin_atom in carboxylic_carbons and mol.GetAtomWithIdx(end_atom).GetAtomicNum() == 8) or
                (end_atom in carboxylic_carbons and mol.GetAtomWithIdx(begin_atom).GetAtomicNum() == 8)
            ):
                continue
            unsaturated_bond = True
            break
    if unsaturated_bond:
        return False, "Unsaturated bonds found in carbon chain"
    
    # Identify carbon chain starting from carboxylic carbon
    from collections import deque

    is_straight_chain = True
    for carboxylic_carbon in carboxylic_carbons:
        visited = set()
        queue = deque()
        queue.append((carboxylic_carbon, -1))  # (current_atom_idx, previous_atom_idx)
        
        while queue:
            current_atom_idx, prev_atom_idx = queue.popleft()
            if current_atom_idx in visited:
                continue
            visited.add(current_atom_idx)
            atom = mol.GetAtomWithIdx(current_atom_idx)
            atomic_num = atom.GetAtomicNum()
            
            # Skip the carboxylic acid oxygens
            if atomic_num != 6 and current_atom_idx not in carboxylic_carbons:
                return False, "Non-carbon atom found in chain (excluding carboxyl group)"
            
            neighbor_idxs = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() != prev_atom_idx]
            if len(neighbor_idxs) > 1:
                # More than one neighbor (excluding previous atom) indicates branching
                if current_atom_idx not in carboxylic_carbons:
                    return False, "Branching detected in carbon chain"
            for neighbor_idx in neighbor_idxs:
                bond = mol.GetBondBetweenAtoms(current_atom_idx, neighbor_idx)
                if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                    return False, "Unsaturated bond detected in carbon chain"
                queue.append((neighbor_idx, current_atom_idx))
    
    return True, "Molecule is a straight-chain saturated fatty acid"

# Example usage (uncomment to test with the provided examples)
# examples = [
#     'CCCCCCCCCCCC(O)=O',  # dodecanoic acid
#     'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O',  # triacontanoic acid
#     'C(CCCCCCCCCCCCCCC(O)=O)CCC(C)O',  # 19-hydroxyicosanoic acid
#     'CCCCCCCCCCCCCCCC(O)=O',  # hexadecanoic acid
#     'OC(C)CCCCCCCCCCCCCCCCCCC(=O)O',  # 20-hydroxyhenicosanoic acid
# ]
# for smiles in examples:
#     result, reason = is_straight_chain_saturated_fatty_acid(smiles)
#     print(f"SMILES: {smiles}\nResult: {result}, Reason: {reason}\n")

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26666',
                              'name': 'straight-chain saturated fatty acid',
                              'definition': 'Any saturated fatty acid lacking a side-chain.',
                              'parents': ['CHEBI:15841']},
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None}