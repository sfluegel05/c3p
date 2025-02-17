"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: a sterol, defined as:
  "Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol 
   (additional carbon atoms may be present in the side chain)."
 
Our approach:
  - Parse the SMILES string.
  - Identify all ring atoms and require that at least 17 of the atoms in the 
    fused ring system are carbons (a rough proxy for the steroid nucleus).
  - Check for the presence of a hydroxyl group (–OH) that is attached to one 
    of the ring carbons.
If both tests are passed, the molecule is regarded as a sterol.
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
    
    # Get ring information; if no rings, cannot be a steroid.
    ring_info = mol.GetRingInfo()
    if ring_info.NumAtomRings() == 0:
        return False, "No rings found in molecule"
    
    # Identify all atoms that participate in any ring.
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
    if not ring_atoms:
        return False, "No ring atoms identified"
    
    # Count the carbon atoms in the ring system.
    carbon_count = 0
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
    if carbon_count < 17:
        return False, ("Fused ring system does not have enough carbon atoms "
                       "to be considered a steroid nucleus (found {}, need >=17)".format(carbon_count))
    
    # Look for hydroxyl groups. We use a SMARTS that matches –OH.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Check if any hydroxyl group is directly attached to a ring carbon.
    found_ring_hydroxyl = False
    for match in hydroxyl_matches:
        o_idx = match[0]  # the oxygen atom of the OH group
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Check neighbors of the hydroxyl oxygen
        for nb in o_atom.GetNeighbors():
            # We require that the hydroxyl is attached to a carbon that is part
            # of the fused ring system.
            if nb.GetAtomicNum() == 6 and nb.GetIdx() in ring_atoms:
                found_ring_hydroxyl = True
                break
        if found_ring_hydroxyl:
            break

    if not found_ring_hydroxyl:
        return False, "No hydroxyl group attached to the ring system found"
    
    return True, ("Molecule has a fused ring system with at least 17 carbon atoms and "
                  "a hydroxyl group attached to one of its ring carbons; "
                  "thus it is classified as a sterol.")

# Example usage:
if __name__ == "__main__":
    # One of the examples provided: (24S,25S)-cholest-5-en-3beta,24,26-triol
    example_smiles = "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_sterol(example_smiles)
    print("Is sterol:", result)
    print("Reason:", reason)