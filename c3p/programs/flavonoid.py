"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid
Definition: A flavonoid is any compound whose skeleton is based on 
a 1-benzopyran core (a fused bicyclic ring system consisting of a heterocycle 
that contains an oxygen fused to an aromatic ring) with an exocyclic aryl substituent 
(i.e. an additional fully aromatic ring attached via a single atom, typically at position 2).

Heuristic implementation:
  1. Parse the molecule using RDKit and obtain all rings.
  2. Look for a pair of fused rings (sharing >=2 atoms) such that one ring contains an oxygen
     (the heterocyclic ring) and at least one of the pair is fully aromatic (the benzene ring).
     Their union is taken as a candidate benzopyran core.
  3. Look for an exocyclic aromatic ring (fully aromatic) that is attached to the candidate core 
     by exactly one atom.
     
Note: Some flavonoids also carry glycoside sugar moieties which should not affect the flavonoid core detection.
This heuristic is not foolproof given the diversity of flavonoid structures.
"""
from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is (putatively) a flavonoid based on its SMILES string.
    
    Heuristic steps:
      1. Parse molecule and get ring data.
      2. Identify a candidate benzopyran core (fused bicyclic system) as a pair of rings
         that are fused (share at least 2 atoms), in which one ring contains an oxygen (to represent the heterocycle)
         and at least one ring is fully aromatic.
      3. Search for an exocyclic aromatic ring (fully aromatic) that attaches to the candidate core
         via exactly one atom (i.e. through a single bond).
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule appears to contain a flavonoid skeleton, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 2:
        return False, f"Only {len(rings)} ring(s) found – expecting at least 2 fused rings for a benzopyran core."

    # Convert rings to sets to simplify intersections.
    ring_sets = [set(r) for r in rings]
    
    candidate_cores = []
    # Loop over all pairs of rings to find one candidate fused pair.
    for i in range(len(ring_sets)):
        for j in range(i+1, len(ring_sets)):
            # We expect fused rings to share at least two atoms.
            if len(ring_sets[i].intersection(ring_sets[j])) < 2:
                continue
            # Check if at least one of the rings is fully aromatic.
            ring1_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_sets[i])
            ring2_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_sets[j])
            # And check that at least one ring contains an oxygen.
            ring1_has_O = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring_sets[i])
            ring2_has_O = any(mol.GetAtomWithIdx(idx).GetSymbol() == 'O' for idx in ring_sets[j])
            if (ring1_has_O or ring2_has_O) and (ring1_aromatic or ring2_aromatic):
                # Create a candidate benzopyran core as the union of the two ring atom indices.
                candidate_core = ring_sets[i].union(ring_sets[j])
                candidate_cores.append(candidate_core)
    
    if not candidate_cores:
        return False, "No fused ring pair found with one ring containing oxygen and one aromatic ring to form a benzopyran core."

    # Now, for each candidate core, search for an exocyclic substituent.
    # The exocyclic aryl substituent must be a fully aromatic ring that is attached to the core by a single atom.
    for core in candidate_cores:
        for ring in ring_sets:
            # Skip rings that are entirely part of the core.
            if ring.issubset(core):
                continue
            # The substituent ring should be fully aromatic.
            if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                continue
            # Determine the intersection between the candidate substituent ring and the core.
            inter = ring.intersection(core)
            if len(inter) == 1:
                return True, "Molecule contains a flavonoid skeleton with a benzopyran core and an exocyclic aryl substituent."
    
    return False, "No appropriate exocyclic fully aromatic substituent attached to a benzopyran core was found."

# Example usage (for testing)
if __name__ == "__main__":
    # Test with (R)-naringenin, a known flavonoid.
    test_smiles = "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    result, reason = is_flavonoid(test_smiles)
    print(result, ":", reason)