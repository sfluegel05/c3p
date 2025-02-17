"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: Saccharolipid – Lipids that contain a carbohydrate moiety.
This function uses two connected heuristics:
  1. It detects a potential sugar ring: here defined as a 5‐ or 6‐membered ring that contains at least one oxygen
     and in which most of the atoms are carbons or substituents typical for sugars (e.g. hydroxyl groups).
  2. For any atom of such a ring, it looks for a “lipid‐like” aliphatic chain attached to it.
     The chain must be a contiguous series of sp³ carbons that are not in any ring. 
     In addition, one bridging oxygen (for instance in an ester or ether linkage) is allowed but only once.
If both conditions are met (with a chain of at least CHAIN_THRESHOLD carbons counted),
the molecule is classified as a saccharolipid.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

# Adjust chain threshold as needed (here we use 10 carbons)
CHAIN_THRESHOLD = 10

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined here as a molecule having a carbohydrate (sugar) ring (ring of size 5 or 6 with at least one O)
    that is connected by a single bond to an aliphatic chain made up of non‐ring sp³ carbons (with at most one bridging oxygen).
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a saccharolipid, False otherwise.
      str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    sugar_found = False
    sugar_atoms = set()
    
    # Look for potential sugar rings: size 5 or 6 with at least one oxygen.
    for ring in ring_info.AtomRings():
        if len(ring) in (5,6):
            oxy_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxy_count += 1
            if oxy_count >= 1:
                sugar_found = True
                sugar_atoms.update(ring)
                
    if not sugar_found:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # For each atom in the sugar ring, look at its neighbors.
    # We only consider neighbors that are not in the sugar ring.
    for idx in sugar_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_atoms:
                continue
            # We only seed DFS from an atom that qualifies: either an aliphatic carbon or a potential bridging oxygen.
            if qualifies_carbon(nbr) or qualifies_oxygen(nbr):
                # Start a DFS that carries a flag 'bridge_used' (False = not yet used)
                chain_len = dfs_chain(nbr, prev_idx=atom.GetIdx(), mol=mol, visited=set(), sugar_set=sugar_atoms, bridge_used=False)
                if chain_len >= CHAIN_THRESHOLD:
                    msg = ("Molecule has a sugar ring with a long aliphatic chain ({} carbon atoms counted) "
                           "attached, classifying it as a saccharolipid".format(chain_len))
                    return True, msg
    return False, "No long aliphatic chain attached to the carbohydrate moiety found"

def qualifies_carbon(atom: rdchem.Atom) -> bool:
    """
    Return True if the atom is a non‐ring sp3 carbon that is not aromatic.
    """
    return (atom.GetAtomicNum() == 6 and
            not atom.GetIsAromatic() and
            not atom.IsInRing() and
            atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)

def qualifies_oxygen(atom: rdchem.Atom) -> bool:
    """
    Return True if the atom is an oxygen (atomic number 8) that is not aromatic and not in a ring.
    (Used as a possible bridge in the chain.)
    """
    return (atom.GetAtomicNum() == 8 and
            not atom.GetIsAromatic() and
            not atom.IsInRing())

def dfs_chain(atom: rdchem.Atom, prev_idx: int, mol: Chem.Mol, visited: set, sugar_set: set, bridge_used: bool) -> int:
    """
    Depth-first search to determine the longest contiguous chain of qualifying aliphatic carbons
    that is reachable from a sugar ring atom. We allow the chain to be “bridged” once through an oxygen.
    
    Args:
      atom: current atom in the DFS.
      prev_idx: index of the atom from which we came (to avoid backtracking).
      mol: the molecule.
      visited: set of already visited atom indices.
      sugar_set: set of indices that are part of sugar rings (do not traverse back into these).
      bridge_used: boolean flag, True if a bridging oxygen has already been used in the current path.
    
    Returns:
      int: maximum count of consecutive carbons (only qualifying carbons count; bridging oxygen contributes 0).
    """
    idx = atom.GetIdx()
    if idx in visited:
        return 0
    if idx in sugar_set:
        return 0  # do not cross back into the sugar moiety
    
    visited.add(idx)
    
    # Determine current contribution: count 1 if a true carbon.
    if qualifies_carbon(atom):
        current = 1
    elif qualifies_oxygen(atom) and (not bridge_used):
        # Allow one bridging oxygen (it contributes no carbons but can link two carbon parts).
        current = 0
        bridge_used = True
    else:
        visited.remove(idx)
        return 0  # stop if the atom does not qualify
    
    max_count = current
    # Look at neighbors (skip the atom we just came from)
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() == prev_idx:
            continue
        if qualifies_carbon(nbr) or qualifies_oxygen(nbr):
            cand = dfs_chain(nbr, prev_idx=idx, mol=mol, visited=visited, sugar_set=sugar_set, bridge_used=bridge_used)
            if current + cand > max_count:
                max_count = current + cand
    visited.remove(idx)
    return max_count

# Example usage:
# Uncomment one of the following test cases to see the classification result.
#
# test_smiles = "[H][C@@]1(O[C@@](C[C@@H](O)[C@H]1O)...)"  # (shortened) example of Lipid A with sugar ring+chain
# print(is_saccharolipid(test_smiles))