"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: Saccharolipid – Lipids that contain a carbohydrate moiety.
This function uses two connected heuristics:
  1. It first detects a sugar ring, defined as a 5‐ or 6‐membered ring with at least one oxygen.
  2. It then looks for a long aliphatic chain. In this improved version the chain is defined as a contiguous
     series of sp³ carbons (counted only if they meet criteria) that can be reached from a sugar ring atom,
     allowing one bridging oxygen (with degree 2) in between.
If both conditions are met the molecule is classified as a saccharolipid.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

# Threshold for the chain (number of aliphatic carbon atoms)
CHAIN_THRESHOLD = 8

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid having a carbohydrate (sugar) moiety
    (a ring of size 5 or 6 with at least one oxygen) that is attached to a long aliphatic chain.
    
    The lipid chain is approximated by detecting a path that yields at least CHAIN_THRESHOLD th sp³ carbon atoms.
    The search is allowed to “bridge” over one oxygen (with degree 2) if it connects two carbon segments.
    
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
    # Step 1. Look for one or more potential sugar rings.
    # A sugar ring is defined here as a 5‐ or 6‐membered ring with at least one oxygen.
    sugar_found = False
    sugar_atoms = set()
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            oxy_count = 0
            for idx in ring:
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
                    oxy_count += 1
            if oxy_count >= 1:
                sugar_found = True
                sugar_atoms.update(ring)
    if not sugar_found:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # Step 2. For each atom in one of the sugar rings, check neighbors
    # to see if a long aliphatic chain (possibly bridged by an oxygen) is attached.
    for idx in sugar_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in sugar_atoms:
                continue  # skip if neighbor is part of the sugar ring
            # Only consider neighbors that are carbon or qualifying oxygen
            if qualifies_carbon(nbr) or qualifies_oxygen(nbr):
                chain_length = dfs_chain(nbr, prev_idx=atom.GetIdx(), mol=mol, visited=set(), sugar_set=sugar_atoms)
                if chain_length >= CHAIN_THRESHOLD:
                    msg = ("Molecule has a sugar ring with a long aliphatic chain (length {} counted as carbons) "
                           "attached, classifying it as a saccharolipid".format(chain_length))
                    return True, msg
    return False, "No long aliphatic chain attached to the carbohydrate moiety found"

def qualifies_carbon(atom: rdchem.Atom) -> bool:
    """
    Returns True if the atom is an aliphatic carbon that is:
      - atomic number 6,
      - not aromatic,
      - not in a ring,
      - sp3-hybridized.
    """
    return (atom.GetAtomicNum() == 6 and
            not atom.GetIsAromatic() and
            not atom.IsInRing() and
            atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)

def qualifies_oxygen(atom: rdchem.Atom) -> bool:
    """
    Returns True if the atom is an oxygen that might serve as a bridge in an ester linkage.
    We require that:
      - the atom is oxygen (atomic number 8),
      - not aromatic,
      - not in a ring,
      - and that it is likely bridging (degree==2).
    """
    return (atom.GetAtomicNum() == 8 and
            not atom.GetIsAromatic() and
            not atom.IsInRing() and
            atom.GetDegree() == 2)

def dfs_chain(atom: rdchem.Atom, prev_idx: int, mol: Chem.Mol, visited: set, sugar_set: set) -> int:
    """
    Using depth-first search we try to determine the maximum count of qualifying aliphatic carbons that form a contiguous chain.
    We permit the chain path to pass through:
      - qualifying carbons (which contribute a count of 1), and
      - qualifying oxygens (which contribute 0 but may bridge parts of the chain).
    The search does not allow revisiting an atom and will not traverse atoms that do not qualify as carbon or bridging oxygen.
    
    Args:
      atom: the current atom in the DFS.
      prev_idx: the index of the atom we came from (to avoid immediate backtracking).
      mol: the molecule.
      visited: a set of already visited atom indices.
      sugar_set: set of atom indices that are part of sugar rings (to avoid counting sugar atoms).
    
    Returns:
      int: the maximum number of carbons that can be counted on a contiguous path from this atom.
    """
    idx = atom.GetIdx()
    if idx in visited:
        return 0
    # Do not cross back into the sugar moiety.
    if idx in sugar_set:
        return 0
    
    visited.add(idx)
    
    # Determine contribution:
    if qualifies_carbon(atom):
        current = 1
    elif qualifies_oxygen(atom):
        current = 0  # oxygen allowed as a bridge
    else:
        visited.remove(idx)
        return 0  # stop the path
    
    max_count = current  # current chain count so far
    
    # Traverse all neighbors (except the one we came from)
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() == prev_idx:
            continue
        # Only allow neighbors that are carbon or a bridging oxygen.
        if qualifies_carbon(nbr) or qualifies_oxygen(nbr):
            cand = dfs_chain(nbr, prev_idx=idx, mol=mol, visited=visited, sugar_set=sugar_set)
            if current + cand > max_count:
                max_count = current + cand
    visited.remove(idx)
    return max_count

# Example usage:
# test_smiles = "CCCCCCCCCCCCCCCC[C@H](C)C[C@H](C)OC1OC(CO)C(O)C(O)C1O"
# print(is_saccharolipid(test_smiles))