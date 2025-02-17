"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: A sulfolipid, defined as 'A compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.'
"""

from rdkit import Chem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is defined as a compound containing a sulfonic acid residue (S(=O)(=O)[O] form)
    that is directly attached (via a C–S bond) to a lipid-like long aliphatic chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a sulfolipid, False otherwise.
        str: Reason for the classification result.
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We'll search for a sulfur atom that is part of a sulfonic acid group
    # i.e. sulfur attached to two oxygens via double bonds and one oxygen via a single bond,
    # and the sulfur must also be attached to a carbon.
    found = False
    reason = "No sulfonic acid residue attached to a lipid chain found"
    
    for atom in mol.GetAtoms():
        # Look for sulfur atoms (atomic number 16)
        if atom.GetAtomicNum() != 16:
            continue
        
        # Get neighboring atoms connected to sulfur
        neighbors = atom.GetNeighbors()
        # To qualify as a sulfonic acid residue, the sulfur should be connected to exactly three atoms
        # (typically 2 oxygens in double bonds and 1 oxygen in a single bond) plus a carbon.
        # However, sometimes the representation can include formal charges or implicit hydrogens,
        # so we will not be overly strict on count but check that at least one neighbor is carbon.
        carbon_neighbor = None
        o_double_count = 0
        o_single_count = 0
        
        for nbr in neighbors:
            # Check for a carbon neighbor
            if nbr.GetAtomicNum() == 6:
                carbon_neighbor = nbr
            # Check oxygen neighbors
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                    o_double_count += 1
                elif bond is not None and bond.GetBondTypeAsDouble() == 1.0:
                    o_single_count += 1
                    
        # We require at least one carbon neighbor and at least two double-bonded oxygens and one single-bonded oxygen.
        # Note: sometimes sulfonate groups are deprotonated and the single-bond oxygen may carry a negative charge.
        if carbon_neighbor is None or o_double_count < 2 or o_single_count < 1:
            continue  # Not the proper sulfonic acid grouping
        
        # Now check whether the carbon neighbor is part of a long aliphatic chain (lipid-like)
        # We do a depth-first search from that carbon to compute the longest contiguous chain length.
        def dfs(atom, visited):
            # Count only carbon atoms that are sp3, non-aromatic and not in a ring
            count = 1
            for nbr in atom.GetNeighbors():
                # Only continue if neighbor is carbon, not already visited,
                # and we restrict to aliphatic (non-aromatic) carbons.
                if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() not in visited):
                    # Check that the neighbor is sp3 and NOT in a ring (typical for lipid chains)
                    if not nbr.GetIsAromatic() and not nbr.IsInRing():
                        new_visited = visited.copy()
                        new_visited.add(nbr.GetIdx())
                        candidate = 1 + dfs(nbr, new_visited)
                        if candidate > count:
                            count = candidate
            return count
        
        # Initialize DFS from the carbon attached to the sulfonic acid S
        chain_length = dfs(carbon_neighbor, {carbon_neighbor.GetIdx()})
        # Use a heuristic cutoff: require at least 8 carbons in a row to qualify as a lipid chain.
        if chain_length >= 8:
            found = True
            reason = f"Found sulfonic acid group attached to a lipid chain (chain length = {chain_length})"
            break
        else:
            reason = f"Found sulfonic acid group, but attached carbon chain too short (chain length = {chain_length})"
    
    if found:
        return True, reason
    else:
        return False, reason