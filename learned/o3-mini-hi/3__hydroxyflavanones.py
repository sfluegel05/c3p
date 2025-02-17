"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3' (meta to the point of attachment on the B ring)
           of the phenyl substituent (B ring). This routine first verifies that a simplified flavanone (chromanone) core is present.
           Then, for each candidate aromatic (six-membered) B ring attached to the core (via the expected C2 atom), we 
           reconstruct ring connectivity to determine the proper meta positions and check that at least one of them carries
           a free hydroxy (-OH) group. If found, the molecule is classified as a member of the class.
           
NOTE: The flavanone core is detected using the SMARTS "C1CC(=O)c2ccccc2O1" which is a simplified version. 
      The ordering of atoms in the B ring is rederived using the graph connectivity (neighbors on the ring) rather than 
      relying on the RDKit ring order.
"""

from rdkit import Chem

def get_ring_neighbors(mol, ring, atom_idx):
    """
    Given a molecule, a set (or list) of atom indices (ring) and a particular atom in that ring,
    return the indices in ring that are directly bonded to the atom.
    """
    neighbors = []
    atom = mol.GetAtomWithIdx(atom_idx)
    for nbr in atom.GetNeighbors():
        if nbr.GetIdx() in ring:
            neighbors.append(nbr.GetIdx())
    return neighbors

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    
    Requirements:
      1. The molecule must contain a flavanone (chromanone) core. For our purpose we use a simplified SMARTS:
         "C1CC(=O)c2ccccc2O1" (this pattern is not perfect but a must-have substructure).
      2. There must be an aromatic six-membered (benzene) ring (B ring) attached to the core.
         We assume that in flavanones the B ring is attached by a single bond to the heterocycle — ideally at the C2.
      3. In that B ring, the atom that is connected to the core is considered as position 1.
         Then its meta positions (those two bonds away along the ring) are determined by re‐computing connectivity.
         At least one of those meta atoms must carry an -OH group (i.e. an oxygen with at least one attached hydrogen)
         that is not part of the benzene ring.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a simplified SMARTS query for the flavanone/chromanone core.
    # (Note: This pattern may capture many flavanone-like cores but is not perfect.)
    core_smarts = "C1CC(=O)c2ccccc2O1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create substructure query for flavanone core"
    
    # Look for the flavanone/chromanone core in the molecule.
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Molecule does not contain a flavanone (chromanone) core"
    
    # We use the first core match found.
    core_match = core_matches[0]
    
    # Look for candidate attachment points: atoms in the core having a neighbor
    # that is both aromatic and NOT in the core.
    candidate_b_attachment = []
    for atom_idx in core_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in core_match and nbr.GetIsAromatic():
                candidate_b_attachment.append(nbr.GetIdx())
    
    if not candidate_b_attachment:
        return False, "Could not find any aromatic substituent (B ring) attached to the flavanone core"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Now, for each candidate attachment, look for a six-membered aromatic ring containing that atom.
    for b_attach in candidate_b_attachment:
        for ring in rings:
            if len(ring) == 6 and b_attach in ring:
                # Check that all atoms in this ring are aromatic (benzene ring).
                if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    continue  # skip if not a proper aromatic ring
                # We now determine the meta positions relative to b_attach.
                # We do not rely on the RDKit ring ordering but recalc connectivity.
                ring_set = set(ring)
                # Get the direct ring neighbors of the attachment atom in the ring.
                direct = get_ring_neighbors(mol, ring_set, b_attach)
                if len(direct) != 2:
                    continue  # something is unusual in this ring; skip
                meta_candidates = set()
                for nbr in direct:
                    nbr_neighbors = get_ring_neighbors(mol, ring_set, nbr)
                    # Exclude the attachment atom itself.
                    for nn in nbr_neighbors:
                        if nn != b_attach:
                            meta_candidates.add(nn)
                if not meta_candidates:
                    continue  # no meta positions found on this ring; try next candidate
                # Check each meta candidate for an -OH substituent.
                for meta_idx in meta_candidates:
                    meta_atom = mol.GetAtomWithIdx(meta_idx)
                    # Look for substituents from this meta atom (neighbors not in the ring).
                    for sub in meta_atom.GetNeighbors():
                        if sub.GetIdx() in ring_set:
                            continue
                        # Check if the substituent is oxygen.
                        if sub.GetAtomicNum() == 8:
                            # Ensure that the oxygen is likely part of a free -OH (has at least one H)
                            if sub.GetTotalNumHs() >= 1:
                                return True, ("Molecule contains flavanone core with a -OH group at the 3' (meta) position " +
                                              "on the B ring")
    # If none of the candidate rings had an appropriate meta –OH substituent.
    return False, "No hydroxy (-OH) group found at the 3' (meta) position on the B ring attached to the flavanone core"

# If this script is run directly, perform a simple test.
if __name__ == "__main__":
    # Test with (2S)-3'-hydroxyflavanone (should be True)
    test_smiles = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, reason = is_3__hydroxyflavanones(test_smiles)
    print(result, reason)