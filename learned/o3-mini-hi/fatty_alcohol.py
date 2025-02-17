"""
Classifies: CHEBI:24026 fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: fatty alcohol
Definition: "An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
             Fatty alcohols may be saturated or unsaturated and may be branched or unbranched."

This implementation:
  1. Checks that the molecule parses.
  2. Rejects molecules with a carboxylic acid/carboxylate group.
  3. Searches for alcohol groups defined as sp³ carbons bound to –OH.
  4. For each candidate alcohol group we require that the alcohol carbon is free – i.e. not in a ring and 
     not attached to any heteroatoms (other than the one –OH).
  5. We then “grow” from that carbon over acyclic, non‐aromatic carbon atoms (i.e. a fatty chain)
     and require that at least one such contiguous chain has a minimum length (here, 7 carbons).
     
If one candidate alcohol passes these checks, we classify the molecule as a fatty alcohol.
"""

from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule belongs to the fatty alcohol class.
    
    The molecule must:
      - Parse correctly.
      - NOT contain a carboxylic acid or carboxylate group.
      - Contain at least one alcohol group (an sp3 carbon bound to –OH) where the alcohol carbon is not in a ring
        and is connected only to carbon neighbors (besides the hydroxyl).
      - Have at least one contiguous, acyclic, non‐aromatic carbon chain (attached to the alcohol carbon) with at 
        least 7 carbons.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple where the boolean is True if the molecule is classified as a fatty alcohol, else False,
                   and the second element is a reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain a carboxylic acid or carboxylate group.
    acid_smarts = "C(=O)[O;H1,-]"  # matches -COOH or -COO-
    acid_query = Chem.MolFromSmarts(acid_smarts)
    if mol.HasSubstructMatch(acid_query):
        return False, "Molecule contains a carboxylic acid/carboxylate group"
    
    # Define SMARTS for an alcohol group: saturated carbon [CX4] attached to hydroxyl (–OH).
    alcohol_smarts = "[CX4][OX2H]"
    alcohol_query = Chem.MolFromSmarts(alcohol_smarts)
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    if not alcohol_matches:
        return False, "No alcohol (-OH) group found on a saturated carbon"
    
    # Build a set of indices for acyclic, non‐aromatic (aliphatic) carbons.
    aliphatic_carbons = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and (not atom.IsInRing()):
            aliphatic_carbons.add(atom.GetIdx())
    if not aliphatic_carbons:
        return False, "No acyclic aliphatic carbon atoms found"
    
    # Build a graph (dictionary) that connects each aliphatic, acyclic carbon to its aliphatic, acyclic neighbors.
    carbon_graph = {idx: [] for idx in aliphatic_carbons}
    for idx in aliphatic_carbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in aliphatic_carbons:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # Helper DFS function to compute the longest simple chain (number of carbons) starting from current node.
    def dfs(current, visited):
        max_length = 1  # count the current atom
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                length = 1 + dfs(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    # Minimum chain length required (e.g., octan-2-ol has 7 carbon chain from the alcohol carbon).
    MIN_CHAIN_LENGTH = 7
    candidate_found = False
    best_chain_length = 0
    
    # Loop over alcohol matches.
    for match in alcohol_matches:
        # match is a tuple, where the first element is the alcohol carbon and the second is the oxygen.
        alc_c_idx = match[0]
        alc_o_idx = match[1]
        alc_atom = mol.GetAtomWithIdx(alc_c_idx)
        # Only consider alcohol carbon if it is not in a ring.
        if alc_atom.IsInRing():
            continue
        
        # Additional filter: the neighbors (other than the hydroxyl oxygen) must be carbon.
        # (This helps exclude cases where the –OH is on a carbon in a decorated environment.)
        neighbors = alc_atom.GetNeighbors()
        valid = True
        for nbr in neighbors:
            # Skip the oxygen that is part of the –OH.
            if nbr.GetIdx() == alc_o_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                valid = False
                break
        if not valid:
            continue

        # Now, check if there is a contiguous chain of aliphatic (acyclic, non‐aromatic) carbons.
        if alc_c_idx not in carbon_graph:
            continue
        chain_length = dfs(alc_c_idx, {alc_c_idx})
        # Track best chain length for reporting.
        if chain_length > best_chain_length:
            best_chain_length = chain_length
        if chain_length >= MIN_CHAIN_LENGTH:
            candidate_found = True
            best_chain_length = chain_length
            # We stop at the first acceptable candidate.
            break

    if not candidate_found:
        return False, f"Longest aliphatic chain starting from an alcohol carbon is only {best_chain_length} atoms (need at least {MIN_CHAIN_LENGTH})"
    
    return True, f"Molecule contains an alcohol group attached to a {best_chain_length}-carbon acyclic aliphatic chain"

# Example usage (for testing):
if __name__ == "__main__":
    # Test with a true positive example: octan-2-ol (SMILES: CCCCCCC(C)O)
    example_smiles = "CCCCCCC(C)O"
    is_fatty, reason = is_fatty_alcohol(example_smiles)
    print(f"SMILES: {example_smiles}\nClassified as fatty alcohol? {is_fatty}\nReason: {reason}")