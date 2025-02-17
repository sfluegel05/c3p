"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: Cholesteryl esters

Definition: A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of cholesterol. In cholesteryl esters the free 3-hydroxy in cholesterol is 
esterified. This improved approach is to first search for an ester-like bond (an O–C(=O) bond) and then, 
for the alcohol-side oxygen, require that it is attached to an aliphatic carbon that is part of a set of 
fused rings whose sizes match the known steroid nucleus (three six-membered rings and one five-membered ring).
"""

from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.

    The method is:
      1. Parse the SMILES string using RDKit.
      2. Iterate over all bonds and select those that are a single bond between an oxygen and a carbon.
      3. Verify that the carbon partner is part of a carbonyl group (i.e. double‐bonded to another oxygen).
      4. For the oxygen (that came from the original 3-hydroxy), check its other neighbor (if any).
         That neighbor should be in a ring. We then collect all rings that contain that neighbor.
      5. We “cluster” rings that share at least two atoms. In cholesterol, the fused steroid nucleus 
         comprises exactly four rings (three of size 6 and one of size 5).
      6. If we find a cluster of exactly four rings whose individual ring sizes (when sorted) equal [5,6,6,6],
         we decide the molecule is a cholesteryl ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): True if the structure is a cholesteryl ester, along with a reason; 
                     False and a reason otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information: list each ring as a set of atom indices.
    ring_info = mol.GetRingInfo()
    all_rings = [set(ring) for ring in ring_info.AtomRings()]
    
    # Helper function: given a list of rings (sets of atom indices),
    # cluster rings that are fused (share at least 2 atoms)
    def cluster_rings(rings):
        clusters = []
        for ring in rings:
            added = False
            for cluster in clusters:
                if any(len(ring & other) >= 2 for other in cluster):
                    cluster.add(frozenset(ring))
                    added = True
                    break
            if not added:
                clusters.append({frozenset(ring)})
        return clusters
    
    # For every bond, search for ester-like O–C(=O) bonds.
    for bond in mol.GetBonds():
        # Require a single bond.
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Identify the pair in which one atom is oxygen and the other is carbon.
        if a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6:
            o_atom = a1
            c_atom = a2
        elif a2.GetAtomicNum() == 8 and a1.GetAtomicNum() == 6:
            o_atom = a2
            c_atom = a1
        else:
            continue
        
        # Check that the carbon (c_atom) is part of a carbonyl group.
        has_carbonyl = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() == o_atom.GetIdx():
                continue  # skip the bond under investigation
            bond_to_nbr = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond_to_nbr.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
        if not has_carbonyl:
            continue  # not an ester bond
        
        # For the alcohol oxygen, check its other neighbor (besides c_atom).
        attached_ring_atom = None
        for nbr in o_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue
            if nbr.IsInRing():
                attached_ring_atom = nbr
                break
        if attached_ring_atom is None:
            continue  # no neighbor in a ring; cannot be a cholesterol derivative.
        
        # Get all rings that include the attached_ring_atom.
        rings_with_atom = [ring for ring in all_rings if attached_ring_atom.GetIdx() in ring]
        if not rings_with_atom:
            continue

        # Cluster these rings according to shared atoms (≥ 2 atoms)    
        clusters = cluster_rings(rings_with_atom)
        # For each fused ring cluster, check if it has the steroid nucleus pattern:
        # exactly 4 rings, with ring sizes 5,6,6,6 (order insensitive)
        for cluster in clusters:
            if len(cluster) != 4:
                continue
            ring_sizes = sorted([len(ring) for ring in cluster])
            if ring_sizes == [5, 6, 6, 6]:
                # Passed the aromaticity and ring-size test typical for cholesterol nucleus.
                return True, ("Contains an ester linkage (O–C(=O)) where the alcohol oxygen is attached to a fused "
                              "steroid nucleus (rings of sizes 5,6,6,6), consistent with a cholesteryl ester")
    return False, ("No ester bond found whose alcohol oxygen is attached to a fused steroid nucleus (four rings with sizes "
                   "5,6,6,6); not a cholesteryl ester")


# Example usage (for testing):
if __name__ == "__main__":
    # Example cholesteryl ester: cholesteryl linoleate
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)"
    result, reason = is_cholesteryl_ester(test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)