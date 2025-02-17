"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3' (meta to the point of attachment
            on the B ring) of the phenyl substituent (B ring). This code first confirms the presence of a 
            simplified flavanone (chromanone) core, then identifies a single six-membered aromatic ring
            (candidate B ring) that attaches to this core at exactly one atom. It then “orders” that ring so
            that the meta positions (two atoms away) are computed and checked for a free –OH.
            
NOTE: This is still a simplified approach. In production further refinement may be required.
"""
from rdkit import Chem

def order_ring(mol, ring_set, start_atom):
    """
    Given a molecule and a set of atom indices that form a ring (assumed to be a simple cycle,
    e.g. benzene), order the ring atoms in sequence starting from start_atom.
    Returns a list of atom indices in order or None if ordering fails.
    """
    ordered = [start_atom]
    current = start_atom
    previous = None
    # For the starting atom, get ring-neighbors within the set:
    neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(current).GetNeighbors() if nbr.GetIdx() in ring_set]
    if not neighbors:
        return None
    # Pick an arbitrary neighbor (our ring is symmetric)
    next_atom = neighbors[0]
    ordered.append(next_atom)
    previous = current
    current = next_atom
    while len(ordered) < len(ring_set):
        nbrs = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(current).GetNeighbors() 
                if nbr.GetIdx() in ring_set and nbr.GetIdx() != previous]
        if not nbrs:
            break
        next_atom = nbrs[0]
        ordered.append(next_atom)
        previous, current = current, next_atom
    if len(ordered) != len(ring_set):
        return None
    return ordered

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    
    Requirements:
      1. The molecule must contain a flavanone (chromanone) core. We use a simplified SMARTS:
         "C1CC(=O)c2ccccc2O1".
      2. The molecule must feature exactly one six-membered aromatic ring (the B ring) that is attached
         to the core by a single bond (a unique attachment atom).
      3. In that benzene ring, when we order the atoms so that the attachment atom is position 0,
         the meta positions (positions 2 and 4 in the ordered list) are examined.
         At least one of these meta carbons must have a free –OH group, defined as an oxygen atom 
         with only one heavy-atom connection (i.e. degree 1 as far as non-hydrogen neighbors) and bearing 
         at least one hydrogen.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule fits the class, False otherwise
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use simplified flavanone (chromanone) core SMARTS.
    core_smarts = "C1CC(=O)c2ccccc2O1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create substructure query for flavanone core"

    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Molecule does not contain a flavanone (chromanone) core"

    # Retrieve ring information for later use.
    ring_info = mol.GetRingInfo().AtomRings()

    # For each flavanone core match, search for the B ring.
    for core_match in core_matches:
        core_set = set(core_match)
        candidate_B_rings = []
        # Review all atoms in the molecule that are aromatic and not part of the core.
        # We want to identify a six-membered aromatic ring that is attached to the core by exactly one bond.
        for ring in ring_info:
            if len(ring) != 6:
                continue
            # Check every atom in this ring is aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            # Count how many atoms in the ring are attached to an atom in the core.
            attach_atoms = []
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in core_set:
                        attach_atoms.append(idx)
                        break
            # We require exactly one atom in the ring is the attachment point.
            if len(attach_atoms) == 1:
                candidate_B_rings.append((ring, attach_atoms[0]))

        if not candidate_B_rings:
            continue  # Try the next core match if none found for this match.

        # For each candidate B ring, check the free -OH condition on the meta positions.
        for ring, attach_idx in candidate_B_rings:
            ring_set = set(ring)
            ordered_ring = order_ring(mol, ring_set, attach_idx)
            if ordered_ring is None:
                continue
            # In an ordered benzene ring, let the attachment index be position 0.
            # The two meta positions (two away along a 6-membered ring) are at indices 2 and 4.
            # (Note: using modulo arithmetic for circularity)
            meta_idx1 = ordered_ring[(ordered_ring.index(attach_idx) + 2) % 6]
            meta_idx2 = ordered_ring[(ordered_ring.index(attach_idx) - 2) % 6]
            meta_positions = [meta_idx1, meta_idx2]
            for meta in meta_positions:
                meta_atom = mol.GetAtomWithIdx(meta)
                # Look for substituents off the meta carbon.
                for nbr in meta_atom.GetNeighbors():
                    # Do not consider atoms in the ring.
                    if nbr.GetIdx() in ring_set:
                        continue
                    # Check if the neighbor is an oxygen.
                    if nbr.GetAtomicNum() == 8:
                        # Check that this oxygen is a free -OH (degree=1 for heavy atoms and at least one H).
                        # (Note: RDKit counts hydrogens automatically unless suppressed.)
                        if nbr.GetDegree() == 1 and nbr.GetTotalNumHs() >= 1:
                            return True, ("Molecule contains flavanone core with a free -OH group at the "
                                          "3' (meta) position on the B ring (attachment at atom idx {} in ring)".format(attach_idx))
        # If none of the candidate B rings satisfy the free -OH condition, try next core match.
    return False, "No free -OH group found at the 3' (meta) position on the B ring attached to the flavanone core"

# When run directly, perform a simple test.
if __name__ == "__main__":
    # Test using (2S)-3'-hydroxyflavanone (should return True)
    test_smiles = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, reason = is_3__hydroxyflavanones(test_smiles)
    print(result, reason)