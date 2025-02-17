"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: Cholesteryl esters

Definition: A sterol ester obtained by formal condensation of a carboxylic acid with the 
3-hydroxy group of cholesterol. In cholesteryl esters the free 3-hydroxy in cholesterol is 
esterified. Here we search for an ester bond in which an oxygen (from the former -OH) is bonded 
to a carbonyl carbon AND that same oxygen is attached to a fused ring system that resembles 
the steroid nucleus (a tetracyclic, fused system). This additional check on “fused‐ness” 
can help weed out false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The method is:
      1. Look for an ester linkage: a single bond between an oxygen and a carbon, where the carbon 
         is double-bonded to another oxygen (i.e. a carbonyl group).
      2. For that ester bond, check that the oxygen (coming from the alcohol part) is attached 
         (via another bond) to a cyclic system that is part of a fused set of rings. In our approach we 
         collect all rings that include that attached atom and then cluster rings that share at least two atoms.
         A fused cluster of 4 or more rings is considered a steroid nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of True with explanation if considered a cholesteryl ester; otherwise False.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get complete ring information for the molecule as sets of atom indices.
    ring_info = mol.GetRingInfo()
    all_rings = [set(ring) for ring in ring_info.AtomRings()]
    
    # Helper function to cluster rings (given as sets) that share at least 2 atoms.
    def cluster_rings(rings):
        clusters = []
        for ring in rings:
            found_cluster = None
            for cluster in clusters:
                # if ring shares at least 2 atoms with any ring in the cluster, add it
                if any(len(ring & other) >= 2 for other in cluster):
                    cluster.add(frozenset(ring))
                    found_cluster = cluster
                    break
            if not found_cluster:
                clusters.append({frozenset(ring)})
        return clusters

    # Iterate over bonds searching for an ester-like linkage
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Identify an oxygen-carbon pair.
        if a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6:
            o_atom = a1
            c_atom = a2
        elif a2.GetAtomicNum() == 8 and a1.GetAtomicNum() == 6:
            o_atom = a2
            c_atom = a1
        else:
            continue
        
        # Check that c_atom forms a carbonyl: it must be double-bonded to an oxygen (other than o_atom)
        has_carbonyl = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() == o_atom.GetIdx():
                continue  # skip the ester-bond oxygen
            bond_to_nbr = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond_to_nbr.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
        if not has_carbonyl:
            continue  # not a proper ester bond
        
        # Now check that the oxygen (o_atom) is also attached to a cyclic system.
        # It will have at least one other neighbor – let’s call that R.
        attached_ring_atom = None
        for nbr in o_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue  # skip the carbonyl part
            # The neighbor must lie in at least one ring.
            if nbr.IsInRing():
                attached_ring_atom = nbr
                break
        if attached_ring_atom is None:
            continue  # the alcohol oxygen is not attached to a cyclic system
        
        # Collect rings that include the attached_ring_atom.
        rings_with_atom = [ring for ring in all_rings if attached_ring_atom.GetIdx() in ring]
        if not rings_with_atom:
            continue  # no ring found for that neighbor (should not occur)
        # Cluster the rings into fused groups (using our simple fusion criteria: share >=2 atoms)
        clusters = cluster_rings(rings_with_atom)
        # Look for any cluster that has at least 4 rings (the steroid nucleus)
        has_steroid_nucleus = any(len(cluster) >= 4 for cluster in clusters)
        if not has_steroid_nucleus:
            continue  # the oxygen is not attached to a sufficiently fused ring system
        
        # If we reach here, we have an ester bond with:
        # - the carbonyl group in place
        # - the oxygen attached (in addition to the carbonyl bond) to an atom that belongs to
        #   a fused ring cluster of at least 4 rings.
        # Thus we consider this as a cholesteryl ester.
        return True, ("Contains an ester linkage (O–C(=O)) with the oxygen attached to a fused cyclic system "
                      "of at least 4 rings (a steroid nucleus), consistent with a cholesteryl ester")
    
    return False, "No ester bond attached to a fused polycyclic (steroid-like) nucleus found; not a cholesteryl ester"


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)"  # cholesteryl linoleate
    result, reason = is_cholesteryl_ester(test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)