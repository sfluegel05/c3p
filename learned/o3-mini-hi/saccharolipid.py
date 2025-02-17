"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: Saccharolipid – Lipids that contain a carbohydrate moiety.
This function uses heuristics:
  1. It finds a sugar ring defined as a 5‐ or 6‐membered ring containing at least one oxygen.
  2. For each atom in such a ring, it searches for a contiguous aliphatic chain (only sp³ carbons)
     attached via a simple bridging oxygen (only one allowed). In addition, the chain must be “clean” –
     meaning that if any atom on the chain is not a simple sp³ carbon (or an allowed bridging oxygen
     that is only connected to carbons), the chain terminates.
If a chain with at least CHAIN_THRESHOLD qualifying carbons is found, the molecule is classified as a saccharolipid.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

# Set the minimum number of aliphatic carbon atoms required in the chain.
CHAIN_THRESHOLD = 10

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as having a carbohydrate (sugar) ring (5‐ or 6‐membered ring with at least one oxygen)
    that is connected (via a single bond) to a long aliphatic chain made up of non‐ring sp³ carbons.
    One bridging oxygen (acting like an ester or ether linkage) is allowed provided the oxygen is only bonded to carbons.
    
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

    # Identify rings of size 5 or 6 that contain at least one oxygen.
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count >= 1:
                sugar_found = True
                sugar_atoms.update(ring)
    
    if not sugar_found:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # For every atom in the sugar ring, check its neighbors for a potential acyl chain.
    for idx in sugar_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # Do not traverse back into the sugar ring.
            if nbr.GetIdx() in sugar_atoms:
                continue
            # Only start from atoms that can be part of an aliphatic chain.
            if qualifies_carbon(nbr) or qualifies_bridge_oxygen(nbr):
                chain_len = dfs_chain(nbr, prev_idx=atom.GetIdx(), mol=mol, sugar_set=sugar_atoms, bridge_used=False, visited=set())
                if chain_len >= CHAIN_THRESHOLD:
                    msg = ("Molecule has a sugar ring with a long aliphatic chain ({} carbon atoms counted) "
                           "attached, classifying it as a saccharolipid".format(chain_len))
                    return True, msg
    return False, "No long aliphatic chain attached to the carbohydrate moiety found"

def qualifies_carbon(atom: rdchem.Atom) -> bool:
    """
    Returns True if the atom is a non‐ring sp3 carbon that is not aromatic.
    """
    return (atom.GetAtomicNum() == 6 and
            not atom.GetIsAromatic() and
            not atom.IsInRing() and
            atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)

def qualifies_bridge_oxygen(atom: rdchem.Atom) -> bool:
    """
    Returns True if the atom is an oxygen that is not aromatic or in a ring
    and appears “clean” (i.e. all its neighbors, except one, are sp3 carbons).
    This function is used to decide if an oxygen can be allowed as a single bridging atom.
    """
    if atom.GetAtomicNum() != 8 or atom.GetIsAromatic() or atom.IsInRing():
        return False

    # Allow oxygen only if (besides the atom we are coming from) all other neighbors are sp3 carbons.
    for nbr in atom.GetNeighbors():
        if nbr.GetAtomicNum() != 6:
            return False
        if nbr.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            return False
    return True

def dfs_chain(atom: rdchem.Atom, prev_idx: int, mol: Chem.Mol, sugar_set: set, bridge_used: bool, visited: set) -> int:
    """
    Depth-first search to determine the maximum number (chain length) of qualifying aliphatic carbons
    along a contiguous chain, starting from an atom attached to the sugar ring.
    Only non‐ring sp3 carbons are counted; a single bridging oxygen (if it qualifies) is allowed
    but contributes 0 to the chain length.
    
    We pass a "visited" set that is copied for each branch so that different DFS branches
    are scored independently.
    
    Args:
      atom: current atom being examined.
      prev_idx: the index of the atom we came from.
      mol: the molecule.
      sugar_set: set of indices belonging to the sugar rings.
      bridge_used: boolean flag indicating whether a bridging oxygen has already been used.
      visited: set of atom indices already visited in the current DFS branch.
    
    Returns:
      int: maximum count of qualifying carbon atoms in a contiguous chain.
    """
    idx = atom.GetIdx()
    if idx in visited:
        return 0
    # Never cross back into the sugar.
    if idx in sugar_set:
        return 0

    # Copy the visited set for current branch (to keep DFS branches independent).
    local_visited = visited.copy()
    local_visited.add(idx)

    current = 0
    if qualifies_carbon(atom):
        current = 1
    elif qualifies_bridge_oxygen(atom) and (not bridge_used):
        # Allow one bridging oxygen (it contributes 0).
        # But before accepting, check that it can connect to at least one qualifying carbon (other than our previous atom).
        valid_next = False
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == prev_idx:
                continue
            if qualifies_carbon(nbr):
                valid_next = True
        if not valid_next:
            return 0
        bridge_used = True
    else:
        return 0  # Atom does not qualify

    max_count = current
    # Explore neighbors (except the one we came from) recursively.
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() == prev_idx:
            continue
        branch_count = dfs_chain(nbr, prev_idx=idx, mol=mol, sugar_set=sugar_set, bridge_used=bridge_used, visited=local_visited)
        if current + branch_count > max_count:
            max_count = current + branch_count
    return max_count

# Example usage:
# Uncomment one of the following test cases to see the classification result.
#
# test_smiles = "[H][C@@]1(O[C@@](C[C@@H](O)[C@H]1O)(O[C@@H]1C[C@@](OC[C@H]2O[C@@H](OC[C@H]3O[C@H](OP(O)(O)=O)..."
# print(is_saccharolipid(test_smiles))