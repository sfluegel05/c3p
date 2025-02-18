"""
Classifies: CHEBI:33566 catechols
"""
#!/usr/bin/env python3
"""
Classifies: Catechols 
Definition: Any compound containing an o-diphenol component.
A catechol is defined as a six‐membered ring (aromatic or highly conjugated)
bearing two adjacent substituents that are hydroxyl groups (–OH or their deprotonated/esterified form),
wherein the oxygen is directly bound to the ring carbon via a single bond and is not further bonded to a carbon.
Moreover, to avoid classifying complex fused polyphenols as “catechols”, we only accept rings that are not fused to any other ring.
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol moiety (an o-diphenol component)
    based on its SMILES string.
    
    The algorithm is as follows:
      1. Convert the SMILES to an RDKit molecule.
      2. Obtain all rings in the molecule. For every six‐membered ring, check if it is “isolated” (not fused) 
         by ensuring that none of its atoms appear in any other ring.
      3. For each isolated ring, first check that the ring is sufficiently conjugated: either the ring is aromatic,
         or – if not aromatic – it contains at least 2 double bonds.
      4. Then, for each bond in that ring (i.e. for every adjacent pair of ring atoms, in an ordered cycle),
         check whether both atoms have a qualifying oxygen substituent. 
         A qualifying oxygen is attached with a single bond and is not further bound to any carbon.
      5. If any pair of adjacent ring atoms meets the criteria, the molecule is classified as a catechol.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a catechol moiety, False otherwise.
        str: Reason for the classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Helper to decide if an oxygen attached to an atom qualifies.
    def has_qualifying_oxygen(atom):
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            # Require bond to be single
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            qualifies = True
            # The oxygen must not be attached to any carbon besides the aromatic/conjugated atom.
            for nb in neighbor.GetNeighbors():
                if nb.GetIdx() == atom.GetIdx():
                    continue
                if nb.GetAtomicNum() == 6:
                    qualifies = False
                    break
            if qualifies:
                return True
        return False

    # Helper to decide whether an atom is present in any other ring than ring_atoms.
    def is_atom_exclusively_in_ring(atom_idx, ring_atoms):
        count = 0
        for ring in atom_rings:
            if atom_idx in ring:
                # if ring is not exactly our candidate ring, count it
                if set(ring) != set(ring_atoms):
                    return False
                count += 1
        return count == 1  # appears only in this ring

    # Helper to check if a candidate ring (list of atom indices) is isolated (non-fused).
    def is_isolated_ring(ring_atoms):
        for idx in ring_atoms:
            # if this atom appears in any ring that is not exactly ring_atoms, consider ring fused.
            if not is_atom_exclusively_in_ring(idx, ring_atoms):
                return False
        return True

    # Helper to order the atoms in a ring.
    # Since atom_rings is an unordered tuple of indices, we produce an ordering by linking adjacent atoms.
    def order_ring(ring_atoms):
        # Build a simple connectivity map restricted to ring_atoms.
        ring_set = set(ring_atoms)
        connectivity = {i: [] for i in ring_atoms}
        for i in ring_atoms:
            atom = mol.GetAtomWithIdx(i)
            for nbr in atom.GetNeighbors():
                j = nbr.GetIdx()
                if j in ring_set:
                    connectivity[i].append(j)
        # Now pick an arbitrary starting atom and perform a simple cycle traversal.
        ordered = [ring_atoms[0]]
        prev = None
        while len(ordered) < len(ring_atoms):
            current = ordered[-1]
            # choose the neighbor that is not the previous atom
            next_atoms = [n for n in connectivity[current] if n != prev]
            if not next_atoms:
                break  # should not happen in a proper cycle
            ordered.append(next_atoms[0])
            prev = current
        return ordered

    # Helper to count the number of explicit double bonds in a ring.
    def count_double_bonds_in_ring(ordered_ring):
        count = 0
        n = len(ordered_ring)
        # because the ring is cyclic, add bond from last to first
        for i in range(n):
            a1 = ordered_ring[i]
            a2 = ordered_ring[(i+1)%n]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                count += 1
        return count

    # Now iterate over all rings of exactly 6 atoms.
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        if not is_isolated_ring(ring):
            # Skip rings that are fused into a larger polycyclic system.
            continue
        ordered_ring = order_ring(list(ring))
        # Determine if the ring is considered aromatic or at least conjugated.
        # If any atom in the ring is flagged as aromatic, treat the ring as aromatic.
        is_aromatic_ring = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ordered_ring)
        if not is_aromatic_ring:
            dbl_bonds = count_double_bonds_in_ring(ordered_ring)
            if dbl_bonds < 2:
                # Not enough conjugation for a catechol moiety.
                continue

        # Now check each adjacent pair in the ring.
        n = len(ordered_ring)
        for i in range(n):
            a1 = mol.GetAtomWithIdx(ordered_ring[i])
            a2 = mol.GetAtomWithIdx(ordered_ring[(i+1)%n])
            # Make sure both atoms are carbons.
            if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
                continue
            # Check both atoms for a qualifying oxygen substituent.
            if has_qualifying_oxygen(a1) and has_qualifying_oxygen(a2):
                return True, "Contains o-diphenol (catechol) moiety on a non-fused, six-membered aromatic or conjugated ring"
                
    return False, "No adjacent qualifying hydroxyl substituents found on an isolated six-membered ring"

# Example usage (for testing):
if __name__ == '__main__':
    examples = [
        # True positives:
        ("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C", "oleuropein"),
        ("C1=CC(=C(C(=C1O)O)[N+]([O-])=O)C", "4-methyl-3-nitrocatechol"),
        ("C1(=CC=C(C=C1O)[N+]([O-])=O)O", "4-nitrocatechol"),
        ("C=1(C=CC(=C(C1)O)O)/C=C/C(OCC)=O", "ethyl trans-caffeate"),
        ("OC[C@H]1O[C@@H](OCCc2ccc(O)c(O)c2)[C@H](OC(=O)Cc2ccc(O)cc2)[C@@H](O)[C@@H]1O", "ternstroside B"),
        ("S(OC1=C(O)C=C([C@@H](O)CN)C=C1)(O)(=O)=O", "Norepinephrine sulfate"),
        ("OC1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1", "flavan-3,3',4,4',5,5',7-heptol"),
        ("COc1cc(-c2ccc(O)cc2)c(OC)c2oc3cc(O)c(O)cc3c12", "candidusin A"),
        ("O=C1OC(=CC2=C1C3(OC(C)=CC3=O)C(C4=CC(O)=C(O)C=C4)O2)/C=C/C5=CC(O)=C(O)C=C5", "Inoscavin A"),
        ("O=C1OC(=CC(=C1C(C=2C(=O)OC(/C=C/C3=CC(O)=C(O)C=C3)=CC2O)C)O)/C=C/C4=CC(O)=C(O)C=C4", "Pinillidine"),
        ("CC(C)c1cccc(O)c1O", "3-isopropylcatechol"),
        ("O[C@@H](CC\\C=C\\c1ccc(O)cc1)CCc1ccc(O)c(O)c1", "(-)-(3S)-1-(3,4-dihydroxyphenyl)-7-(4-hydroxyphenyl)-(6E)-6-hepten-3-ol"),
        # False negative example (should be detected):
        ("O[C@H]1C=CC=C[C@H]1O", "cis-cyclohexa-3,5-diene-1,2-diol"),
    ]
    for smi, name in examples:
        result, reason = is_catechols(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {result} -> {reason}\n{'-'*80}")