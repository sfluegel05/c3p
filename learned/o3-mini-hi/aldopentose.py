"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: aldopentose, defined as 'A pentose with a (potential) aldehyde group at one end.'
An aldopentose may appear as an open‐chain sugar – having a terminal aldehyde – 
or as a cyclic (hemiacetal) sugar in equilibrium with the open‐chain form.
We no longer simply count all C and O atoms, since cyclic forms like pyranoses 
may give misleading totals. Instead we extract the sugar “backbone” by finding 
the longest connected chain of carbon atoms (ignoring extra carbons from substituents).
For an aldopentose this longest chain should have 5 carbons.
Furthermore, if an aldehyde functionality (SMARTS "[CX3H1](=O)") is present, we check 
that it is on a terminal carbon of that backbone.
If not (but the molecule is cyclic and does not contain a lactone/ester "[CX3](=O)[O]"),
we classify it as a cyclized aldopentose.
"""

from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    
    An aldopentose sugar (open‐chain form C5H10O5) should ultimately have 5 carbons in its sugar 
    backbone. In the open‐chain form a terminal carbon should bear an aldehyde group.
    In cyclic (furanose or pyranose) forms, although the total atom counts can be misleading,
    the longest contiguous chain of carbons (ignoring oxygen in the ring) is still 5.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an aldopentose; False otherwise.
        str: Explanation for the classification decision.
    """    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Build a subgraph containing only carbon atoms.
    # We create a mapping from carbon atom index to its neighboring carbon indices.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return False, "No carbon atoms found in molecule"
    carbon_graph = {idx: set() for idx in carbon_idxs}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].add(idx2)
                carbon_graph[idx2].add(idx1)
    
    # Use DFS to compute the maximum length of a simple path within the carbon subgraph.
    # (The length is counted as the number of carbon atoms in that path.)
    def dfs(current, visited):
        max_len = 1
        for neighbor in carbon_graph[current]:
            if neighbor not in visited:
                new_len = 1 + dfs(neighbor, visited | {neighbor})
                if new_len > max_len:
                    max_len = new_len
        return max_len

    longest_chain = 0
    for idx in carbon_idxs:
        chain_length = dfs(idx, {idx})
        if chain_length > longest_chain:
            longest_chain = chain_length

    if longest_chain != 5:
        return False, f"Longest carbon chain length is {longest_chain}; expected 5 for a pentose sugar backbone"
    
    # Look for an aldehyde group.
    # The SMARTS "[CX3H1](=O)" will match an sp2 carbon attached to one hydrogen and double-bonded to O.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    
    # For a proper open‐chain aldopentose, the aldehyde should reside at a terminal carbon on the sugar backbone.
    open_chain_aldehyde = False
    # We use the carbon_graph built above: a terminal carbon in that subgraph will have degree 1.
    for match in aldehyde_matches:
        # match[0] is the carbon atom in the aldehyde group.
        carbon_idx = match[0]
        # Check if this carbon is in our carbon backbone graph and if it is terminal.
        if carbon_idx in carbon_graph and len(carbon_graph[carbon_idx]) == 1:
            open_chain_aldehyde = True
            break

    # Also define lactone (ester) pattern so we can exclude lactone sugars.
    lactone_pattern = Chem.MolFromSmarts("[CX3](=O)[O]")
    
    if open_chain_aldehyde:
        return True, "Open‐chain aldopentose: terminal aldehyde group detected on the sugar backbone."
    else:
        # If no open-chain aldehyde, then for a cyclic sugar we expect a ring (furanose or pyranose)
        # and no lactone (ester) functionality.
        if mol.GetRingInfo().NumRings() > 0:
            if mol.HasSubstructMatch(lactone_pattern):
                return False, "Cyclic structure contains lactone functionality; not an aldopentose."
            return True, "Cyclized aldopentose: no free aldehyde but cyclic, in equilibrium with open‐chain form."
        else:
            return False, "Molecule is acyclic and lacks a terminal aldehyde group; does not meet aldopentose criteria."
    
# (Optional: you may include testing calls below this line.)