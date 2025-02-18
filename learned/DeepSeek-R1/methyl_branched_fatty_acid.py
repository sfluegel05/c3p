"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
"""
Classifies: methyl-branched fatty acid (CHEBI: 84878)
A methyl-branched fatty acid is a branched-chain fatty acid containing methyl branches only.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, MolFromSmarts

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    Must have a carboxylic acid group and only methyl branches on the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid = MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Get the main chain starting from carboxylic acid
    # Find longest carbon chain starting from COOH alpha carbon
    cooh_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not cooh_matches:
        return False, "No carboxylic acid group"
    
    # Alpha carbon is the one attached to COOH
    alpha_carbons = {match[1] for match in cooh_matches}  # [O]=C-O connects to alpha carbon at index 1
    
    # Function to find longest chain from alpha carbon
    def find_longest_chain(start_atom_idx):
        visited = set()
        stack = [(start_atom_idx, 0)]  # (atom_idx, chain_length)
        max_length = 0
        while stack:
            atom_idx, length = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            max_length = max(max_length, length + 1)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in alpha_carbons:
                    stack.append((neighbor.GetIdx(), length + 1))
        return max_length

    # Find the longest possible chain from each alpha carbon
    max_chain_length = max(find_longest_chain(a) for a in alpha_carbons)
    
    # Now check all carbons in the main chain for branches
    # For simplicity, assume main chain is correctly identified as longest possible
    # Check all carbons not in COOH group for non-methyl branches
    for atom in mol.GetAtoms():
        if atom.GetIdx() in alpha_carbons:
            continue  # skip alpha carbon (part of COOH)
        if atom.GetAtomicNum() != 6:
            continue
        
        # Check for branching: if atom has more than 2 carbons attached (excluding main chain)
        neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(neighbors) > 2:
            return False, f"Branching at carbon {atom.GetIdx()+1} with degree >2"
        
        # Check substituents (non-main chain carbons)
        for neighbor in neighbors:
            # Check if neighbor is part of the main chain
            # This part is simplified; may need more accurate chain tracking
            # Instead, check if substituent has more than one carbon (i.e., not methyl)
            substituent = neighbor
            sub_carbons = 0
            stack = [substituent]
            visited = set()
            while stack:
                current = stack.pop()
                if current.GetIdx() in visited:
                    continue
                visited.add(current.GetIdx())
                if current.GetAtomicNum() == 6:
                    sub_carbons += 1
                    if sub_carbons > 1:
                        return False, f"Non-methyl substituent at carbon {atom.GetIdx()+1}"
                    for n in current.GetNeighbors():
                        if n.GetIdx() != atom.GetIdx():
                            stack.append(n)
    
    # Verify at least one methyl branch exists
    methyl_pattern = MolFromSmarts("[CH3]")
    has_methyl = False
    for match in mol.GetSubstructMatches(methyl_pattern):
        methyl_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in methyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in alpha_carbons:
                has_methyl = True
                break
        if has_methyl:
            break
    if not has_methyl:
        return False, "No methyl branches found"

    return True, "Contains only methyl branches on fatty acid chain"