"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 4'-hydroxyflavanones
Definition: “Any hydroxyflavanone having a hydroxy substituent located at position 4'.”
Heuristic: 
  (1) The molecule must contain a cyclic ketone substructure (to indicate the flavanone chromanone core).
  (2) It must contain an isolated aromatic six‐membered ring (the B ring) that is attached via one of its
      atoms to a non‐aromatic fragment (i.e. the flavanone core).
  (3) In that aromatic ring, one of the atoms should bear an –OH substituent in the para position relative 
      to the attachment point (i.e. separated by three bonds along the ring cycle).
"""

from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if the molecule represented by `smiles` is a 4'-hydroxyflavanone. 
    It first checks for the presence of a cyclic ketone (as a proxy for the flavanone core),
    then iterates over aromatic six-membered rings (candidate B rings) to see if one
    has an -OH group attached para (i.e. 1,4 relative) to the atom linking it to the core.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if molecule is classified as a 4'-hydroxyflavanone, False otherwise.
       str: A reason explaining the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for cyclic ketone substructure.
    # We use a SMARTS to roughly find a carbonyl (C=O) that is part of a ring.
    ketone_smarts = Chem.MolFromSmarts("[R]C(=O)[R]")
    if not mol.HasSubstructMatch(ketone_smarts):
        return False, "No cyclic ketone (flavanone core) detected"

    # 2. Get all ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each is a tuple of atom indices
    
    # For each aromatic six-membered ring (candidate for the flavanone B-ring):
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # we want six-membered aromatic rings
        # Verify that every atom in the ring is aromatic
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # Try to identify the attachment point.
        # Look for an atom in this ring that is connected to an atom NOT in the ring.
        attachment_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors not in the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring:
                    # To be safe, we consider this an attachment if the neighbor is not aromatic
                    # (typically the flavanone core is not aromatic at the attachment point)
                    if not nbr.GetIsAromatic():
                        attachment_indices.append(idx)
                        break
        if not attachment_indices:
            # No attachment found on this ring; likely not the B-ring
            continue

        # Now, look for an -OH group substituent on one of the ring atoms.
        # For each atom in the ring, check if it is bound to an oxygen that is part of an -OH unit.
        hydroxy_positions = []  # will store indices in the ring (positions defined by the order in 'ring')
        # We assume the order in the 'ring' tuple gives a cyclic order.
        for pos, idx in enumerate(ring):
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Check that the neighbor is an oxygen and is not part of the ring.
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                    # Check that this oxygen is an -OH (it should have at least one hydrogen)
                    # Note: H atoms are often implicit but RDKit can report implicit H count.
                    if nbr.GetTotalNumHs() > 0:
                        hydroxy_positions.append(pos)
                        break  # found an -OH on this ring atom; move to next ring atom

        if not hydroxy_positions:
            continue  # no -OH groups on this ring, try next ring

        # For each identified attachment atom (its position within the ring tuple)...
        for attach_idx in attachment_indices:
            try:
                pos_attach = ring.index(attach_idx)
            except ValueError:
                continue  # should not happen
            # Now check if any -OH group in the ring is located para relative to the attachment.
            for pos_oh in hydroxy_positions:
                # In a six-membered ring, the para position relative to a given atom is 3 bonds away.
                # Since the ring is a cycle, use the minimal distance between positions along the cycle.
                dist = abs(pos_oh - pos_attach)
                if dist > 3:  # adjust for wrap-around 
                    dist = 6 - dist
                if dist == 3:
                    return True, ("Molecule has a cyclic ketone core and a six-membered aromatic (B) ring "
                                  "with an -OH substituent para to the attachment (position 4').")
                    
    return False, "No appropriate six-membered aromatic B ring with a para -OH (4'-OH) substituent detected"

# (Optional) Testing examples:
if __name__ == '__main__':
    # A few examples are given in the prompt (e.g. taxifolin)
    test_smiles = [
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin (should be True)
        "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",             # sakuranetin (should be True)
        "Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1",                  # naringenin (should be True)
        "COC(=O)C1=CC=CC=C1"                                   # An ester, not a flavanone (should be False)
    ]
    for sm in test_smiles:
        result, reason = is_4__hydroxyflavanones(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")