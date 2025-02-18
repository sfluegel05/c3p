"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: Neoflavonoid 
Definition: A neoflavonoid is any 1-benzopyran (chromene) that has an aryl (e.g. phenyl) substituent at its 4-position.
The core is defined as a fused bicyclic system of a 6-membered heterocycle (chromene, containing one oxygen)
fused to a benzene ring. The aryl substituent must be attached at a carbon that is not part of the fused (shared) region.
"""

from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid.
    A neoflavonoid must have:
      1. A 6-membered heterocyclic ring (chromene) containing exactly one oxygen.
      2. A benzene ring (6 aromatic carbons) fused (sharing at least 2 atoms) with the heterocycle.
      3. An aryl substituent attached to one of the non-fused positions of the heterocycle (heuristically the 4-position).
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is a neoflavonoid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    chromene_ring = None  # candidate 6-membered heterocycle containing exactly one oxygen
    benzene_ring = None   # candidate benzene ring: 6 atoms, all aromatic carbons
    
    # Search for candidate rings
    for ring in rings:
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            # Check for candidate benzene: all aromatic and atomic number=6 (C)
            if all(a.GetIsAromatic() and a.GetAtomicNum() == 6 for a in atoms):
                benzene_ring = ring
            # Check for candidate heterocycle (chromene part): exactly one O and at least one C.
            if sum(1 for a in atoms if a.GetAtomicNum() == 8) == 1 and sum(1 for a in atoms if a.GetAtomicNum() == 6) >= 1:
                chromene_ring = ring
    
    if chromene_ring is None:
        return False, "No 6-membered heterocycle (potential chromene) containing one oxygen was found"
    if benzene_ring is None:
        return False, "No benzene ring (aromatic 6-membered carbon ring) was found"
    
    # Check that the heterocycle and benzene ring are fused (share at least 2 atoms)
    shared_atoms = set(chromene_ring).intersection(set(benzene_ring))
    if len(shared_atoms) < 2:
        return False, "Candidate chromene and benzene rings are not properly fused"
    
    # In the chromene ring, determine non-fused atoms (i.e. not shared with the benzene ring)
    non_fused_atoms = set(chromene_ring) - shared_atoms
    if not non_fused_atoms:
        return False, "No non-fused positions found on the heterocycle"
    
    # Look for an external aryl substituent on any non-fused carbon atom
    # (We interpret this as an external neighbor that is aromatic and belongs to a 6-membered aromatic ring.)
    for idx in non_fused_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # Consider only carbon atoms (position 4 should be a carbon)
        if atom.GetAtomicNum() != 6:
            continue
        # Iterate over neighbors not in the chromene ring
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in chromene_ring:
                continue  # skip atoms that are part of the core
            # Check if this neighbor is aromatic
            if not nbr.GetIsAromatic():
                continue
            # Now check if the neighbor is part of a 6-membered aromatic ring (e.g. phenyl ring)
            for ring in ring_info.AtomRings():
                if nbr.GetIdx() in ring and len(ring) == 6:
                    # Verify that every atom in the ring is aromatic
                    if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                        return True, "Found chromene core fused with benzene and an aryl substituent at a candidate non-fused position"
                        
    return False, "No suitable aryl substituent attached at a non-fused position on the chromene core could be identified"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with one of the provided example SMILES strings.
    test_smiles = "CC1=CC(=CC2=C1C(CC3(O2)CC(NC(=S)N3)(C)C)C4=CC=C(C=C4)OC)O"
    result, reason = is_neoflavonoid(test_smiles)
    print("Is neoflavonoid?", result)
    print("Reason:", reason)