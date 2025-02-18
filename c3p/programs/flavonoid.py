"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid 
Definition: Any member of the 'superclass' flavonoids whose skeleton is based on 
1-benzopyran with an aryl substituent at position 2.

This implementation does not simply search for a 2-phenylbenzopyran SMARTS. Instead, 
it uses RDKit’s ring information to first require that the molecule contains at least 
three rings, then identifies (a) aromatic rings (candidate flavonoid ring “A” and “B”) 
and (b) candidate heterocyclic rings (flavonoid “C” containing an oxygen). Next it requires 
that one heterocycle is fused (sharing at least two atoms) with an aromatic ring and that 
it bears an exocyclic aromatic substituent. This exocyclic substituent is taken 
as consistent with the aryl group at position 2.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a (putative) flavonoid based on its SMILES string.
    The heuristic first checks that the molecule has at least three rings.
    Then it identifies a candidate heterocyclic ring (ring C) that:
       - contains an oxygen and is not fully aromatic,
       - is fused (shares 2 or more atoms) with an aromatic ring (ring A) 
         (thus forming the classical benzopyran core),
       - and that has at least one exocyclic aromatic neighbor (ring B) attached
         (consistent with an aryl substituent at position 2).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to have a flavonoid skeleton; False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # tuple of tuples (each tuple is a ring atom indices)
    
    # Flavonoids normally have three rings.
    if len(rings) < 3:
        return False, f"Only {len(rings)} ring(s) found – expecting at least 3 for a flavonoid skeleton."
    
    # Collect aromatic rings (which are fully aromatic)
    aromatic_rings = []
    for ring in rings:
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(set(ring))
    
    # Collect candidate heterocyclic rings: not fully aromatic but containing at least one oxygen.
    candidate_heterocycles = []
    for ring in rings:
        ring_set = set(ring)
        has_oxygen = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring)
        is_fully_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if has_oxygen and (not is_fully_aromatic):
            candidate_heterocycles.append(ring_set)
    
    if not candidate_heterocycles:
        return False, "No candidate heterocycle (containing oxygen and not fully aromatic) found."
    
    # Now try to find a candidate heterocycle that is fused to an aromatic ring.
    # Fused means sharing at least two atoms.
    for hetero_ring in candidate_heterocycles:
        for aro_ring in aromatic_rings:
            if len(hetero_ring.intersection(aro_ring)) >= 2:
                # Now, from any atom in the candidate heterocycle, see if we have an exocyclic aromatic substituent.
                # Exocyclic means a neighbor that is not part of the heterocycle.
                for atom_idx in hetero_ring:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for nbr in atom.GetNeighbors():
                        nbr_idx = nbr.GetIdx()
                        if nbr_idx not in hetero_ring:
                            # Check that this neighbor is part of some aromatic ring (i.e. a candidate for ring B).
                            # We require the neighbor itself to be aromatic.
                            if nbr.GetIsAromatic():
                                return True, "Molecule contains a flavonoid skeleton with a fused benzopyran core and an aryl substituent at position 2."
    
    return False, "Molecule does not appear to contain a flavonoid skeleton as defined."

# Example usage (for testing)
if __name__ == "__main__":
    # Test with (R)-naringenin: a widely known flavonoid.
    test_smiles = "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    result, reason = is_flavonoid(test_smiles)
    print(result, ":", reason)