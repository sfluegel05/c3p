"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3' (meta to the point of attachment
            on the B ring) of the phenyl substituent (B ring). This version first confirms the presence of a
            flavanone (chromanone) core (using a simplified SMARTS pattern), then identifies the unique aromatic
            substituent attached (ideally from the C2 position). It then “orders” that benzene ring so that the meta
            positions (two atoms away from the connection) are exactly identified, and finally verifies that one
            of those meta carbons carries a free –OH group (an oxygen attached to that carbon only and to hydrogens).
            
NOTE: Even with these improvements, using a simplified core SMARTS ("C1CC(=O)c2ccccc2O1") may pick up extra
      structures. Further refinement might be needed for production use.
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
    # There must be two ring-neighbors, pick one arbitrarily as next:
    neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(current).GetNeighbors() if nbr.GetIdx() in ring_set]
    if not neighbors:
        return None
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
      2. The B ring (an aromatic six-membered ring) must be attached to the core at the expected carbon (ideally the C2 atom).
      3. In that benzene ring, if we designate the attachment atom as position 1 then the meta positions
         (positions 3 and 5 in the ordered benzene cycle) are computed and at least one must have a free –OH group.
         A free –OH group is defined as an oxygen atom that is attached only to that ring carbon (degree 1 on heavy atoms)
         and has at least one hydrogen.
    
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
    
    # Create substructure query for simplified flavanone/chromanone core.
    core_smarts = "C1CC(=O)c2ccccc2O1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create substructure query for flavanone core"
    
    # Find matches of the core in the molecule.
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Molecule does not contain a flavanone (chromanone) core"
    
    # For each match in the core, we try to locate the attached aromatic (B ring).
    # We assume the attachment should come from a non-aromatic (i.e. aliphatic) atom in the core.
    for core_match in core_matches:
        # For each atom in the core match that is non-aromatic:
        for atom_idx in core_match:
            core_atom = mol.GetAtomWithIdx(atom_idx)
            if core_atom.GetIsAromatic():
                continue  # We want the saturated part of the chromanone core.
            # Look for neighbors not in the core that are aromatic. This is a candidate for the B ring attachment.
            for nbr in core_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in core_match or not nbr.GetIsAromatic():
                    continue
                # Now, verify that this candidate lies in a six-membered aromatic ring.
                rings = mol.GetRingInfo().AtomRings()
                candidate_rings = []
                for ring in rings:
                    if len(ring) == 6 and nbr_idx in ring:
                        # Ensure every atom in the ring is aromatic.
                        if all(mol.GetAtomWithIdx(r).GetIsAromatic() for r in ring):
                            candidate_rings.append(ring)
                if not candidate_rings:
                    continue
                # For each valid candidate ring, try to order the ring atoms.
                for ring in candidate_rings:
                    ring_set = set(ring)
                    ordered_ring = order_ring(mol, ring_set, nbr_idx)
                    if ordered_ring is None:
                        continue
                    # In the ordered benzene, define the attachment position as index 0.
                    # The meta positions in a benzene ring (6-membered) are at indices 2 and 4.
                    meta_indices = [ordered_ring[(ordered_ring.index(nbr_idx) + 2) % 6],
                                    ordered_ring[(ordered_ring.index(nbr_idx) - 2) % 6]]
                    # Check if either meta carbon has a free -OH group.
                    for meta_idx in meta_indices:
                        meta_atom = mol.GetAtomWithIdx(meta_idx)
                        # For each neighbor off the meta carbon that is not in the benzene ring:
                        for sub in meta_atom.GetNeighbors():
                            if sub.GetIdx() in ring_set:
                                continue
                            if sub.GetAtomicNum() == 8:
                                # Confirm that the oxygen is free:
                                # It should be attached only to this meta carbon (degree 1 in heavy atoms)
                                # and have at least one hydrogen.
                                if sub.GetDegree() == 1 and sub.GetTotalNumHs() >= 1:
                                    return True, ("Molecule contains flavanone core with a free -OH group at the "
                                                  "3' (meta) position on the B ring")
    return False, ("No free hydroxy (-OH) group found at the 3' (meta) position on the B ring attached to the flavanone core")

# When run directly, perform a simple test.
if __name__ == "__main__":
    # Test using (2S)-3'-hydroxyflavanone (should return True)
    test_smiles = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, reason = is_3__hydroxyflavanones(test_smiles)
    print(result, reason)