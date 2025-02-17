"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: Saccharolipid – Lipids that contain a carbohydrate moiety.
This improved function uses two connected heuristics:
  1. It first detects a sugar ring, defined as a 5‐ or 6‐membered ring with at least one oxygen.
  2. It then looks for a long aliphatic chain (here at least 8 connected sp3 carbons that are
     non‐aromatic, non‐ring) that is attached to one of the sugar ring atoms.
If both conditions are met (and the chain is “attached” to the sugar), the molecule is classified as a saccharolipid.
"""

from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid having a carbohydrate (sugar) moiety (a ring of size 5 or 6
    with at least one oxygen) that is attached to a long aliphatic chain.
    
    The lipid part is approximated by detecting a chain of at least 8 contiguous sp³ carbons that
    are non‐aromatic, not in a ring, and directly attached to a sugar atom.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a saccharolipid, False otherwise.
      str: Reason for the classification.
    """
    # parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()

    # First, search for one or more potential sugar rings.
    # We require rings of size 5 or 6 where at least one atom is oxygen.
    sugar_found = False
    sugar_atoms = set()  # will hold indices of atoms in the recognized sugar rings
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            oxygen_count = 0
            for idx in ring:
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
                    oxygen_count += 1
            if oxygen_count >= 1:
                sugar_found = True
                sugar_atoms.update(ring)
    if not sugar_found:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # Next, for each atom in the sugar ring(s), check its neighbors to see
    # if a long aliphatic (lipid) chain is attached.
    # The chain is defined as a contiguous series of sp3-hybridized carbons (atomic number 6)
    # that are non-aromatic and not in any ring (they are “free” chains).
    for idx in sugar_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # We require the neighbor to be a carbon, non-aromatic, and not already in a ring.
            if (nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() and not nbr.IsInRing()):
                # Start a DFS from the neighbor to find the longest contiguous aliphatic chain length.
                chain_length = _dfs_longest_chain(nbr, mol, {nbr.GetIdx()}, sugar_atoms)
                if chain_length >= 8:
                    return True, ("Molecule has a sugar ring with a long aliphatic chain (length {} "
                                  "attached), classifying it as a saccharolipid".format(chain_length))
    
    return False, "No long aliphatic chain attached to the carbohydrate moiety found"

def _dfs_longest_chain(atom, mol, visited, sugar_set):
    """
    Depth-first search to compute the length of the longest chain from a starting atom.
    Only traverses carbon atoms that are: 
      - sp3-hybridized (typical for aliphatic chains),
      - non-aromatic,
      - not in a ring,
      - and not part of any recognized sugar (sugar_set).
    Parameters:
      atom (rdkit.Chem.rdchem.Atom): starting atom.
      mol (rdkit.Chem.Mol): molecule.
      visited (set): set of already visited atom indices.
      sugar_set (set): indices that belong to the sugar ring(s) (to avoid counting those).
    Returns:
      int: the chain length (number of atoms in the longest contiguous path).
    """
    max_length = 1  # count the current atom
    for nbr in atom.GetNeighbors():
        nbr_idx = nbr.GetIdx()
        # only continue if neighbor hasn't been visited, is a carbon, non-aromatic, not in a ring,
        # and not part of the sugar moiety.
        if (nbr_idx not in visited and 
            nbr.GetAtomicNum() == 6 and 
            not nbr.GetIsAromatic() and 
            not nbr.IsInRing() and
            nbr_idx not in sugar_set):
            visited.add(nbr_idx)
            length = 1 + _dfs_longest_chain(nbr, mol, visited, sugar_set)
            if length > max_length:
                max_length = length
            visited.remove(nbr_idx)
    return max_length

# Example usage:
# test_smiles = "CCCCCCCCCCCCCCCC[C@H](C)C[C@H](C)OC1OC(CO)C(O)C(O)C1O" 
# print(is_saccharolipid(test_smiles))