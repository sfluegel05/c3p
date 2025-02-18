"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 3'-hydroxyflavanones
Definition: Any hydroxyflavanone with a hydroxy substituent at position 3' (meta to the point of attachment on the B
            ring) of the phenyl substituent (B ring). This routine first verifies that a simplified flavanone (chromanone)
            core is present. Then, for each candidate aromatic (six-membered) B ring attached to the core (via the expected
            C2 atom), we reconstruct ring connectivity to determine the proper meta positions and check that at least one
            of them carries a free hydroxy (-OH) substituent. The additional check now verifies that the oxygen in the -OH
            group is only attached to the benzene ring (i.e. it is not further substituted as in a glycoside or ester).
            
NOTE: The flavanone core is detected using the SMARTS "C1CC(=O)c2ccccc2O1" which is a simplified version.
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
      1. The molecule must contain a flavanone (chromanone) core. Here we use a simplified SMARTS:
         "C1CC(=O)c2ccccc2O1". (This pattern may match extra structures so caution is needed.)
      2. There must be an aromatic six-membered (benzene) ring (B ring) attached to the core.
      3. In that B ring, the atom that is connected to the core is considered as position 1.
         Its meta positions (two bonds away in the ring) are recomputed using connectivity.
         At least one of those meta atoms must carry a free –OH group (an oxygen that is attached only
         to that ring carbon and to a hydrogen), not involved in a glycosidic/ester link.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise
        str: Reason for classification.
    """
    # Parse the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Create substructure query for the simplified flavanone/chromanone core.
    core_smarts = "C1CC(=O)c2ccccc2O1"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create substructure query for flavanone core"
    
    # Look for the flavanone/chromanone core in the molecule.
    core_matches = mol.GetSubstructMatches(core_query)
    if not core_matches:
        return False, "Molecule does not contain a flavanone (chromanone) core"
    
    # We use the first core match.
    core_match = core_matches[0]
    
    # Look among atoms in the core for a neighbor that is aromatic and outside the core.
    candidate_b_attachment = []
    for atom_idx in core_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in core_match and nbr.GetIsAromatic():
                candidate_b_attachment.append(nbr.GetIdx())
    
    if not candidate_b_attachment:
        return False, "Could not find any aromatic substituent (B ring) attached to the flavanone core"
    
    # Get ring information (all rings) in the molecule.
    rings = mol.GetRingInfo().AtomRings()
    
    # For each candidate attachment, check if they lie in a six-membered aromatic ring.
    for b_attach in candidate_b_attachment:
        for ring in rings:
            if len(ring) != 6 or b_attach not in ring:
                continue
            # Verify that all atoms in this ring are aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            # Define the ring as a set for faster membership checks.
            ring_set = set(ring)
            # Get the immediate ring neighbors of the attachment atom.
            direct_neighbors = get_ring_neighbors(mol, ring_set, b_attach)
            if len(direct_neighbors) != 2:
                continue  # Unexpected ring connectivity; skip this ring.
            meta_candidates = set()
            # For each direct neighbor, its other neighbor (besides the attachment) is considered meta.
            for nbr in direct_neighbors:
                nbr_neighbors = get_ring_neighbors(mol, ring_set, nbr)
                for nn in nbr_neighbors:
                    if nn != b_attach:
                        meta_candidates.add(nn)
            if not meta_candidates:
                continue  # No meta positions found in this ring.
            # Check each meta candidate for a substituent that is a free –OH group.
            for meta_idx in meta_candidates:
                meta_atom = mol.GetAtomWithIdx(meta_idx)
                # Look at the substituents off the meta atom, ignoring atoms inside the ring.
                for sub in meta_atom.GetNeighbors():
                    if sub.GetIdx() in ring_set:
                        continue
                    # Check if the substituent is oxygen.
                    if sub.GetAtomicNum() == 8:
                        # Ensure that the oxygen is a free –OH.
                        # A free –OH oxygen should be attached only to the meta atom (heavy atom count == 1)
                        # and have at least one hydrogen.
                        if sub.GetDegree() == 1 and sub.GetTotalNumHs() >= 1:
                            return True, ("Molecule contains flavanone core with a free -OH group at the 3' (meta) position "
                                          "on the B ring")
    return False, "No free hydroxy (-OH) group found at the 3' (meta) position on the B ring attached to the flavanone core"

# When run directly, perform a simple test.
if __name__ == "__main__":
    # Test using (2S)-3'-hydroxyflavanone (expected True)
    test_smiles = "Oc1cccc(c1)[C@@H]1CC(=O)c2ccccc2O1"
    result, reason = is_3__hydroxyflavanones(test_smiles)
    print(result, reason)