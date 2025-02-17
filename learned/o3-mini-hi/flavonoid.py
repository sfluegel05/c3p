"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid
Definition: A flavonoid is defined as any compound whose skeleton is based on 
a 1-benzopyran core (i.e. a fused aromatic+heterocyclic ring system where the 
heterocycle contains at least one oxygen) with an aryl substituent (an exocyclic 
aromatic ring) at position 2.
    
This implementation uses RDKit’s ring information. It first checks that the 
molecule has at least three rings. Then it collects (a) aromatic rings (candidate “ring B” 
for the aryl substituent) and (b) oxygen-containing rings (candidate “ring C”). 
For each oxygen-containing ring that is fused (shares at least two atoms) with an aromatic ring 
(thereby forming a potential benzopyran core), the code then checks if any atom on the core 
has an exocyclic aromatic neighbor. The exocyclic aromatic neighbor is defined as belonging to 
an aromatic ring that shares exactly one atom with the fused core.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a (putative) flavonoid based on its SMILES string.
    
    Heuristic:
      - The molecule must have at least three rings.
      - One ring (the candidate heterocycle) must contain at least one oxygen.
      - That oxygen ring must be fused (share at least two atoms) with an aromatic ring 
        (thus forming the benzopyran core).
      - The fused core must have at least one exocyclic aromatic substituent 
        (an aromatic ring attached through a single atom with that aromatic ring otherwise 
         disjoint from the core).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to contain a flavonoid skeleton; False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    if len(rings) < 3:
        return False, f"Only {len(rings)} ring(s) found – expecting at least 3 for a flavonoid skeleton."
    
    # Identify all fully aromatic rings (candidate for the exocyclic aryl substituents)
    aromatic_rings = []
    for ring in rings:
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    
    # Identify rings that contain at least one oxygen atom (candidate for the heterocycle part of the benzopyran)
    oxygen_rings = []
    for ring in rings:
        ring_set = set(ring)
        if any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring):
            oxygen_rings.append(ring_set)
    
    if not oxygen_rings:
        return False, "No ring containing oxygen found; cannot form a benzopyran core."
    
    # Now look for a candidate oxygen ring (ring C) that is fused with an aromatic ring.
    for ox_ring in oxygen_rings:
        for aro_ring in aromatic_rings:
            # Fused means sharing at least 2 atoms.
            if len(ox_ring.intersection(aro_ring)) >= 2:
                # Define the candidate benzopyran core as the union of the oxygen ring and the fused aromatic ring.
                candidate_core = ox_ring.union(aro_ring)
                # Look for an exocyclic aromatic substituent attached to the candidate core.
                # Exocyclic: a neighbor of a core atom that is part of an aromatic ring that 
                # does not merge with the core (its intersection with the core is a single atom).
                for atom_idx in candidate_core:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for nbr in atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx not in candidate_core:
                            # Check if this neighbor is in an aromatic ring that is connected via only one atom.
                            for ext_ring in aromatic_rings:
                                if nbr_idx in ext_ring:
                                    # If the ring is attached to the core by exactly one shared atom,
                                    # we treat it as the aryl substituent.
                                    if len(ext_ring.intersection(candidate_core)) == 1:
                                        return True, ("Molecule contains a flavonoid skeleton with a fused benzopyran core "
                                                      "and an aryl substituent at position 2.")
    return False, "Molecule does not appear to contain a flavonoid skeleton as defined."

# Example usage (for testing)
if __name__ == "__main__":
    # Example test SMILES for (R)-naringenin, a well-known flavonoid.
    test_smiles = "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    result, reason = is_flavonoid(test_smiles)
    print(result, ":", reason)