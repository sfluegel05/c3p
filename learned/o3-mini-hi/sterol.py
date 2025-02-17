"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: a sterol, defined as:
  "Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol 
   (additional carbon atoms may be present in the side chain)."

Our approach:
  - Parse the SMILES string.
  - Verify that the molecule contains rings by checking the AtomRings list.
  - Determine the fused ring system by taking the union of all ring atom indices.
  - Count the carbon atoms in the ring system; require at least 17 (proxy for steroid nucleus).
  - Check for a hydroxyl group (-OH) attached directly to one of the ring carbons.
"""

from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol.
    A sterol is defined as any 3-hydroxy steroid having a tetracyclic ring system
    (closely related to cholestan-3-ol) with a hydroxyl group attached to a ring carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a sterol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the ring information. Instead of using NumAtomRings(), we use AtomRings().
    atom_rings = mol.GetRingInfo().AtomRings()
    if not atom_rings or len(atom_rings) == 0:
        return False, "No rings found in molecule"
    
    # Identify all atoms that participate in any ring (the fused ring system).
    ring_atoms = set()
    for ring in atom_rings:
        ring_atoms.update(ring)
    if not ring_atoms:
        return False, "No ring atoms identified"
    
    # Count the number of carbon atoms in the ring system.
    carbon_count = 0
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
    # This threshold is a rough cutoff for a steroid nucleus.
    if carbon_count < 17:
        return False, ("Fused ring system does not have enough carbon atoms "
                       f"to be considered a steroid nucleus (found {carbon_count}, need >= 17)")
    
    # Look for hydroxyl groups; use a SMARTS that matches -OH (hydroxyl group).
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Check if any hydroxyl group is directly attached to a ring carbon.
    found_ring_hydroxyl = False
    for match in hydroxyl_matches:
        # match is a tuple with the oxygen atom index from the hydroxyl group.
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Check neighbors of the hydroxyl oxygen.
        for nb in o_atom.GetNeighbors():
            # If the neighboring atom is carbon and is in the ring system, we count it.
            if nb.GetAtomicNum() == 6 and nb.GetIdx() in ring_atoms:
                found_ring_hydroxyl = True
                break
        if found_ring_hydroxyl:
            break

    if not found_ring_hydroxyl:
        return False, "No hydroxyl group attached to the fused ring system was found"
    
    # Passed both conditions: fused steroid nucleus and hydroxyl on ring carbon.
    return True, ("Molecule has a fused ring system with at least 17 carbon atoms and "
                  "a hydroxyl group attached to one of its ring carbons; "
                  "thus, it is classified as a sterol.")

# Example usage:
if __name__ == "__main__":
    # Test one of the example SMILES for a sterol
    test_smiles = "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_sterol(test_smiles)
    print("Is sterol:", result)
    print("Reason:", reason)