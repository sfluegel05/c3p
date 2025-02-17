"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: Cholesteryl esters

Definition: A sterol ester obtained by formal condensation of the carboxy group of any carboxylic acid 
with the 3-hydroxy group of cholesterol. In cholesteryl esters the free 3-hydroxy in cholesterol is 
esterified. Our approach is to first search for an ester linkage (an O–C(=O) bond) and then examine 
the oxygen (which came originally from cholesterol) to see if it is attached to a fused cyclic system 
that is of a size typical for the cholesterol nucleus (roughly 15–20 carbon atoms in the fused rings).
"""

from rdkit import Chem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    
    The method is:
      1. Parse the SMILES to get an RDKit molecule.
      2. Iterate over all bonds; if a bond is a single bond between an oxygen and a carbon,
         check that the carbon (other than the oxygen) is double‐bonded to another oxygen (i.e. part
         of a carbonyl group).
      3. For that ester bond, examine the oxygen (the “alcohol oxygen”) – it should have 
         one more neighbor. If that neighbor is part of a ring, collect all rings containing that atom, 
         and cluster them using a simple criterion: rings that share at least 2 atoms.
      4. If a fused ring cluster is found, count the number of carbon atoms in the union; if this number 
         is in the range 15 to 20 (typical for the cholesterol nucleus) we decide that the compound is a 
         cholesteryl ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple where True indicates a cholesteryl ester with a reason; False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information as sets of atom indices.
    ring_info = mol.GetRingInfo()
    all_rings = [set(ring) for ring in ring_info.AtomRings()]
    
    # Helper function: cluster rings that share at least two atoms.
    def cluster_rings(rings):
        clusters = []
        for ring in rings:
            found_cluster = None
            for cluster in clusters:
                if any(len(ring & other) >= 2 for other in cluster):
                    cluster.add(frozenset(ring))
                    found_cluster = cluster
                    break
            if not found_cluster:
                clusters.append({frozenset(ring)})
        return clusters

    # Iterate over bonds looking for an ester-like linkage.
    for bond in mol.GetBonds():
        # Require a single bond (the ester O-C bond).
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        
        # Identify a pair in which one atom is oxygen and the other is carbon.
        if a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6:
            o_atom = a1
            c_atom = a2
        elif a2.GetAtomicNum() == 8 and a1.GetAtomicNum() == 6:
            o_atom = a2
            c_atom = a1
        else:
            continue
        
        # Check that the carbon (c_atom) is carbonyl: it must be double-bonded to an oxygen
        has_carbonyl = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() == o_atom.GetIdx():
                continue  # skip the bond just under consideration
            bond_to_nbr = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8 and bond_to_nbr.GetBondType() == Chem.BondType.DOUBLE:
                has_carbonyl = True
                break
        if not has_carbonyl:
            continue  # not an ester bond
        
        # Now check the oxygen (o_atom) for its second attachment.
        # It should be attached (beside c_atom) to an atom in a ring.
        attached_ring_atom = None
        for nbr in o_atom.GetNeighbors():
            if nbr.GetIdx() == c_atom.GetIdx():
                continue
            if nbr.IsInRing():
                attached_ring_atom = nbr
                break
        if attached_ring_atom is None:
            continue  # The ester oxygen is not attached to any cyclic system.
        
        # Collect all rings that include the attached_ring_atom.
        rings_with_atom = [ring for ring in all_rings if attached_ring_atom.GetIdx() in ring]
        if not rings_with_atom:
            continue
        
        clusters = cluster_rings(rings_with_atom)
        # In the fused steroid nucleus we expect a cluster of 4 (or so) rings.
        # But we further check that the total number of carbon atoms in the unified set of atoms 
        # in one of those clusters falls into a typical range (15 to 20) for the cholesteryl nucleus.
        for cluster in clusters:
            # Merge all atom indices in the cluster.
            cluster_atoms = set()
            for ring in cluster:
                cluster_atoms |= set(ring)
            # Count carbon atoms in that cluster.
            carbon_count = sum(1 for idx in cluster_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if 15 <= carbon_count <= 20:
                return True, ("Contains an ester linkage (O–C(=O)) where the oxygen is attached to a fused cyclic system "
                              f"with {carbon_count} carbons, consistent with a cholesteryl ester")
    return False, ("No appropriate ester bond attached to a fused polycyclic nucleus (with 15–20 carbons) found; "
                   "not a cholesteryl ester")


# Example usage (for testing):
if __name__ == "__main__":
    # Example cholesteryl ester: cholesteryl linoleate
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@H](CC[C@]4(C)[C@@]3([H])CC[C@]12C)OC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)"
    result, reason = is_cholesteryl_ester(test_smiles)
    print("Classification result:", result)
    print("Reason:", reason)