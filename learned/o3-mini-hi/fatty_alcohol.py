"""
Classifies: CHEBI:24026 fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: fatty alcohol
Definition: "An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
            Fatty alcohols may be saturated or unsaturated and may be branched or unbranched."
            
This function attempts to classify a molecule as a fatty alcohol by:
  1. Parsing the SMILES.
  2. Checking for at least one alcohol group ([CX4][OX2H]) – i.e. an sp³ carbon attached to –OH.
  3. Building a connectivity graph of aliphatic (non‐aromatic) carbon atoms.
  4. For any alcohol carbon in an alcohol group, computing the longest connected chain (by DFS) 
     among aliphatic carbons.
  5. If a chain of at least three carbons is found, we return True with the appropriate reason.
  
If the analysis cannot be done the function will return (False, <reason>).
"""

from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) belongs to the fatty alcohol class.
    
    For our purposes, the molecule must:
      - Be valid and parseable.
      - Contain at least one alcohol group (an sp3 carbon bound to a hydroxyl group).
      - Possess a contiguous aliphatic (non‐aromatic) carbon chain (connected via C–C bonds)
        of at least 3 atoms that is connected to an alcohol carbon.
        
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a fatty alcohol, False otherwise.
      str: A reason explaining the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS for an alcohol group: a sp3 (saturated) carbon attached to a hydroxyl.
    alcohol_smarts = "[CX4][OX2H]"
    alcohol_query = Chem.MolFromSmarts(alcohol_smarts)
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    if not alcohol_matches:
        return False, "No alcohol (–OH) group found on a saturated carbon"

    # Build a list of indices for aliphatic (non-aromatic) carbon atoms.
    aliphatic_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.GetIsAromatic():
            aliphatic_carbons.append(atom.GetIdx())

    if not aliphatic_carbons:
        return False, "No aliphatic carbon atoms found"

    # Build a dictionary mapping carbon atom index -> list of neighboring carbon indices (only aliphatic carbons).
    carbon_graph = {idx: [] for idx in aliphatic_carbons}
    for idx in aliphatic_carbons:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() in aliphatic_carbons):
                carbon_graph[idx].append(nbr.GetIdx())

    # Helper function: DFS to compute the longest simple path starting from a given carbon.
    def dfs(current, visited):
        max_length = 1  # count the current atom
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                length = 1 + dfs(nbr, visited | {nbr})
                if length > max_length:
                    max_length = length
        return max_length

    # For each alcohol match, check whether the carbon (the first atom in the match) is in a long chain.
    # We require a chain length of at least 3 (per definition). Note that the chain may be branched,
    # but we are looking for the longest connected sub-chain in the aliphatic region.
    chain_found = False
    max_chain_length_found = 0
    for match in alcohol_matches:
        alcohol_c_idx = match[0]  # index of the carbon attached to the OH.
        if alcohol_c_idx not in carbon_graph:
            continue
        # DFS starting from alcohol carbon.
        chain_length = dfs(alcohol_c_idx, {alcohol_c_idx})
        if chain_length > max_chain_length_found:
            max_chain_length_found = chain_length
        if chain_length >= 3:
            chain_found = True
            # We stop early if one alcohol group is connected to a sufficient chain.
            break

    if not chain_found:
        return False, f"Longest aliphatic chain starting from an alcohol carbon is only {max_chain_length_found} atoms (need at least 3)"
    
    # Optionally: one could include further checks such as checking the overall molecular size or purity
    # of the aliphatic chain. Our current criteria are based only on the presence of an alcohol group
    # attached to a sufficiently long aliphatic chain.
    return True, f"Molecule contains an alcohol group attached to a {max_chain_length_found}-carbon aliphatic chain"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples: octan-2-ol: SMILES: CCCCCCC(C)O
    example_smiles = "CCCCCCC(C)O"
    is_fatty, reason = is_fatty_alcohol(example_smiles)
    print(f"SMILES: {example_smiles}\nClassified as fatty alcohol? {is_fatty}\nReason: {reason}")