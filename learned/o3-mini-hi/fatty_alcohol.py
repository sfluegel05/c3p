"""
Classifies: CHEBI:24026 fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: fatty alcohol
Definition: "An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
             Fatty alcohols may be saturated or unsaturated and may be branched or unbranched."
             
This improved version uses several additional criteria:
  1. The molecule must be parseable.
  2. It must NOT contain a carboxylic acid/carboxylate group.
  3. It must contain at least one alcohol group defined by an sp³ carbon (i.e. [CX4])–OH.
  4. The alcohol carbon must be “exposed” (i.e. not inside a ring) and be attached to
     a contiguous chain of aliphatic (non‐aromatic, non‐ring) carbons having a minimum length.
  5. We require a connected chain of at least 7 carbons (this minimum helps filter out small‐sized alcohols).
     
If one such alcohol group is found with an acceptable chain length, the molecule is classified as a fatty alcohol.
"""

from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) belongs to the fatty alcohol class.
    
    The molecule must:
      - Be valid and parseable.
      - NOT contain a carboxylic acid group.
      - Contain at least one alcohol group (an sp3 carbon bound to an –OH) that is not in a ring.
      - Have a contiguous, aliphatic (non‐aromatic, acyclic) carbon chain attached to the alcohol carbon 
        that is at least 7 carbon atoms long.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a fatty alcohol, False otherwise.
      str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject if molecule contains a carboxylic acid or carboxylate group.
    # SMARTS matches both acid (–COOH) and carboxylate (–COO–).
    acid_smarts = "C(=O)[O;H1,-]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    if mol.HasSubstructMatch(acid_query):
        return False, "Molecule contains a carboxylic acid/carboxylate group"

    # Define SMARTS for an alcohol group: sp3 (saturated) carbon attached to –OH:
    alcohol_smarts = "[CX4][OX2H]"
    alcohol_query = Chem.MolFromSmarts(alcohol_smarts)
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    if not alcohol_matches:
        return False, "No alcohol (–OH) group found on a saturated carbon"
    
    # Build a list of indices for aliphatic (non‐aromatic) carbons that are also acyclic.
    # We will use these to build connectivity among carbons that very likely belong to fatty chains.
    aliphatic_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic() and not atom.IsInRing():
            aliphatic_carbons.append(atom.GetIdx())
    if not aliphatic_carbons:
        return False, "No acyclic aliphatic carbon atoms found"
        
    # Build a graph: for each aliphatic, acyclic carbon, list its neighboring carbons (also meeting same criteria).
    carbon_graph = {idx: [] for idx in aliphatic_carbons}
    for idx in aliphatic_carbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() in aliphatic_carbons):
                carbon_graph[idx].append(nbr.GetIdx())
    
    # Helper DFS function to compute the longest simple path (in terms of number of carbons) starting from a given node.
    def dfs(current, visited):
        max_length = 1   # count the current atom
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                length = 1 + dfs(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    # Set a minimum carbon chain length for a fatty alcohol. (For example, octan-2-ol should have a chain of 7.)
    MIN_CHAIN_LENGTH = 7

    chain_found = False
    max_chain_length_found = 0
    # Loop over each alcohol match.
    for match in alcohol_matches:
        alcohol_c_idx = match[0]  # the alcohol group is defined as [CX4][OX2H]; take the carbon part.
        # Only consider alcohol carbons that are acyclic (in our aliphatic graph).
        atom = mol.GetAtomWithIdx(alcohol_c_idx)
        if atom.IsInRing():
            continue  # skip alcohols whose carbons lie in a ring (likely not part of a free fatty chain).
        if alcohol_c_idx not in carbon_graph:
            continue
        chain_length = dfs(alcohol_c_idx, {alcohol_c_idx})
        if chain_length > max_chain_length_found:
            max_chain_length_found = chain_length
        if chain_length >= MIN_CHAIN_LENGTH:
            chain_found = True
            # Stop early if one acceptable alcohol is found.
            break
    
    if not chain_found:
        return False, f"Longest aliphatic chain starting from an alcohol carbon is only {max_chain_length_found} atoms (need at least {MIN_CHAIN_LENGTH})"

    return True, f"Molecule contains an alcohol group attached to a {max_chain_length_found}-carbon acyclic aliphatic chain"

# Example usage:
if __name__ == "__main__":
    # Test with one true positive example: octan-2-ol (SMILES: CCCCCCC(C)O)
    example_smiles = "CCCCCCC(C)O"
    is_fatty, reason = is_fatty_alcohol(example_smiles)
    print(f"SMILES: {example_smiles}\nClassified as fatty alcohol? {is_fatty}\nReason: {reason}")