"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid
Definition: A flavonoid is defined as any compound whose skeleton is based on 
a 1-benzopyran core (i.e. a fused bicyclic ring system consisting of a heterocycle 
that contains an oxygen fused to an aromatic ring) with an aryl substituent 
(i.e. a fully aromatic ring attached via a single atom) at position 2.
    
This implementation uses RDKit’s ring information. It first requires that the 
molecule has at least three rings (the benzopyran core plus the exocyclic aryl ring). 
Then it searches for any two rings that are fused (i.e. share at least two atoms) 
where one ring contains at least one oxygen and the other is fully aromatic 
(allowing for flavanones where the heterocycle is not fully aromatic). 
If such a fused pair is found (candidate benzopyran core) and an additional fully aromatic 
ring is attached to the core (by sharing exactly one atom), the molecule is classified as a flavonoid.
    
Note: This is a heuristic – flavonoid structures are quite diverse and some edge cases may be missed.
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a (putative) flavonoid based on its SMILES string.
    
    Heuristic steps:
      1. Parse the molecule and obtain ring information.
      2. Require that there are at least 3 rings (the fused benzopyran core and the exocyclic aryl ring).
      3. Look for a fused pair of rings (sharing at least 2 atoms) where one ring
         contains oxygen and the other ring is fully aromatic (this pair forms the flavonoid benzopyran core).
      4. Then, check for an additional ring that is fully aromatic and attached to the fused core
         via exactly one atom (this is the aryl substituent at position 2).
         
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule appears to contain a flavonoid skeleton, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information: each ring is returned as a tuple of atom indices.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    if len(rings) < 3:
        return False, f"Only {len(rings)} ring(s) found – expecting at least 3 for a flavonoid skeleton."
    
    # Convert rings to sets for easier set operations.
    ring_sets = [set(r) for r in rings]
    
    candidate_cores = []
    # Look for any two rings that are fused (i.e. share at least two atoms).
    # For the flavonoid benzopyran core, we require one of the rings to contain oxygen
    # and the other ring to be fully aromatic.
    for i in range(len(ring_sets)):
        for j in range(i+1, len(ring_sets)):
            intersection = ring_sets[i].intersection(ring_sets[j])
            if len(intersection) < 2:
                continue  # not fused, skip
            # Check if one ring contains oxygen and the other is fully aromatic.
            ring1_has_O = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring_sets[i])
            ring2_has_O = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring_sets[j])
            ring1_is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_sets[i])
            ring2_is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_sets[j])
            
            if (ring1_has_O and ring2_is_aromatic) or (ring2_has_O and ring1_is_aromatic):
                # form candidate benzopyran core from the union of these two rings
                candidate_core = ring_sets[i].union(ring_sets[j])
                candidate_cores.append(candidate_core)
    
    if not candidate_cores:
        return False, "No fused ring pair with an oxygen-containing ring and a fully aromatic ring detected to form a benzopyran core."
    
    # Now search for an exocyclic aromatic substituent:
    # An exocyclic ring should be fully aromatic and attach to the candidate core via exactly one atom.
    for core in candidate_cores:
        for ring in ring_sets:
            # Skip if the ring is already part of the core.
            if ring.issubset(core):
                continue
            # We require the exocyclic substituent to be fully aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            # Check if the intersection with the core is exactly one atom.
            if len(ring.intersection(core)) == 1:
                return True, "Molecule contains a flavonoid skeleton with a benzopyran core and an aryl substituent at position 2."
    
    return False, "No appropriate exocyclic aryl substituent attached to a benzopyran core was found."

# Example usage (for testing)
if __name__ == "__main__":
    # Test SMILES for (R)-naringenin, a known flavonoid.
    test_smiles = "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    result, reason = is_flavonoid(test_smiles)
    print(result, ":", reason)