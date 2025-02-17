"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Cyclic Fatty Acid
Definition: Any fatty acid containing anywhere in its structure a ring of atoms.
A fatty acid here is expected to have a carboxylic acid moiety, a long (>=5 carbons) contiguous acyclic
aliphatic chain, and at least one ring present somewhere in the structure.
"""

from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    
    A cyclic fatty acid is defined as having:
      - a carboxylic acid group (matches both protonated and deprotonated forms),
      - at least one ring anywhere in its structure, and 
      - a contiguous chain of 5 or more non-ring, non-aromatic carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a cyclic fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for carboxylic acid.
    # This pattern covers both the protonated ([OX2H1]) and deprotonated ([OX1-]) forms.
    carboxylic_acid_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H1,OX1-]')
    if carboxylic_acid_pattern is None:
        return False, "Failed to generate carboxylic acid SMARTS pattern"
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check if the molecule contains at least one ring.
    ring_info = mol.GetRingInfo()
    if not ring_info or len(ring_info.AtomRings()) == 0:
        return False, "No ring found in the structure"

    # Instead of using a SMARTS with a repeating operator, we now compute the longest chain of 
    # non-ring, non-aromatic (i.e., aliphatic) carbon atoms.
    # Identify indices of such carbons.
    acyclic_carbons = [atom.GetIdx() for atom in mol.GetAtoms() 
                         if atom.GetAtomicNum() == 6 and not atom.IsInRing() and not atom.GetIsAromatic()]
                         
    # Build adjacency for these acyclic carbons.
    neighbors = {idx: [] for idx in acyclic_carbons}
    for idx in acyclic_carbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in acyclic_carbons:
                neighbors[idx].append(nbr.GetIdx())
    
    # Use DFS to compute the longest contiguous chain in the acyclic carbon subgraph.
    def dfs(current, visited):
        max_length = 1
        for nbr in neighbors[current]:
            if nbr not in visited:
                length = 1 + dfs(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    longest_chain = 0
    for idx in acyclic_carbons:
        chain_length = dfs(idx, {idx})
        if chain_length > longest_chain:
            longest_chain = chain_length
    
    if longest_chain < 5:
        return False, "No sufficiently long contiguous acyclic aliphatic chain (>=5 carbons) found"
    
    return True, f"Contains carboxylic acid group, a contiguous acyclic carbon chain of length {longest_chain}, and a ring in the structure"

# For debugging and testing:
if __name__ == "__main__":
    # Example: (R)-lipoic acid is a known cyclic fatty acid.
    test_smiles = "OC(=O)CCCC[C@@H]1CCSS1"
    result, reason = is_cyclic_fatty_acid(test_smiles)
    print(result, reason)