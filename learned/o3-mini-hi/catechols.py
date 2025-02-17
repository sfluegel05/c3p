"""
Classifies: CHEBI:33566 catechols
"""
#!/usr/bin/env python3
"""
Classifies: Catechols
Definition: Any compound containing an o-diphenol component.
A catechol moiety is defined here as a six‐membered ring (aromatic or sufficiently conjugated)
that bears two adjacent, qualifying hydroxyl (-OH) substituents.
A “qualifying” hydroxyl is one that is bound to an aromatic (sp2 carbon) via a single bond
and is not further linked to another carbon (i.e. not part of an ether, ester, etc.).
Note:
  – We do not enforce “isolation” (i.e. non-fusion with any other ring) because that rule was found to miss cases.
  – This classifier therefore attempts to capture many of the examples where an o-diphenol substructure exists.
  
There remain trade‐offs since many natural products contain catechol sub‐units embedded in fused or poly‐hydroxylated systems.
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule contains a catechol moiety (an o-diphenol component)
    based on its SMILES string.
    
    The algorithm is as follows:
      1. Convert the SMILES to an RDKit molecule.
      2. Retrieve all rings; for every six‐membered ring, ensure that it is aromatic or, if not,
         contains at least 2 explicit double bonds.
      3. For each candidate ring, order the ring atoms and, for each aromatic carbon, check if it bears
         a “qualifying” hydroxyl substituent. A qualifying OH is one attached by a single bond
         and not further bound to any carbon.
      4. If any pair of adjacent carbons (in the ring order) both carry a qualifying –OH group,
         we classify the molecule as containing a catechol moiety.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a catechol moiety, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to get a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices

    # Helper function: decide if a given oxygen, bound to an atom in the ring,
    # qualifies as a hydroxyl (–OH) rather than (for example) an ether substituent.
    def has_qualifying_oh(atom):
        # Look at each neighbor of the atom.
        # We require an oxygen (atomic number 8) attached via a single bond.
        # Then, check that the oxygen is not attached to any carbon except for “atom”.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Count the number of neighboring carbons (other than the current atom).
            # If none, then we assume this oxygen is present as an isolated –OH.
            attached_carbons = sum(1 for x in nbr.GetNeighbors() if x.GetAtomicNum() == 6 and x.GetIdx() != atom.GetIdx())
            if attached_carbons == 0:
                return True
        return False

    # Helper function: return an ordering of the ring atoms.
    # Since ring_info returns an unordered tuple, we “walk” the ring using connectivity.
    def order_ring(ring_atoms):
        ring_set = set(ring_atoms)
        connectivity = {i: [] for i in ring_atoms}
        for i in ring_atoms:
            atom = mol.GetAtomWithIdx(i)
            for nbr in atom.GetNeighbors():
                j = nbr.GetIdx()
                if j in ring_set:
                    connectivity[i].append(j)
        ordered = [ring_atoms[0]]
        prev = None
        while len(ordered) < len(ring_atoms):
            current = ordered[-1]
            # choose neighbor different from the previous atom
            next_atoms = [n for n in connectivity[current] if n != prev]
            if not next_atoms:
                break  # in theory should not happen if ring is connected
            ordered.append(next_atoms[0])
            prev = current
        return ordered

    # Helper function: decide if a ring (given by its ordered list of atom indices) is conjugated.
    # We consider it conjugated if any atom is flagged aromatic,
    # or if not aromatic, if the ring contains at least 2 double bonds.
    def ring_is_conjugated(ordered_ring):
        atoms = [mol.GetAtomWithIdx(i) for i in ordered_ring]
        if any(atom.GetIsAromatic() for atom in atoms):
            return True
        dbl = 0
        n = len(ordered_ring)
        for i in range(n):
            a1 = ordered_ring[i]
            a2 = ordered_ring[(i+1)%n]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                dbl += 1
        return dbl >= 2

    # Now consider each six-membered ring.
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        ordered_ring = order_ring(list(ring))
        if not ring_is_conjugated(ordered_ring):
            continue

        # For each atom in the ring, record if it is an aromatic carbon with a qualifying –OH.
        qualifier = {}
        for idx in ordered_ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider aromatic carbons (for our catechol pattern)
            if atom.GetAtomicNum() == 6 and atom.GetIsAromatic():
                qualifier[idx] = has_qualifying_oh(atom)
            else:
                qualifier[idx] = False

        # Check for an adjacent pair in the ordered ring that both have qualifying –OH groups.
        n = len(ordered_ring)
        for i in range(n):
            if qualifier[ordered_ring[i]] and qualifier[ordered_ring[(i+1)%n]]:
                # If found, we report success.
                return True, "Contains o-diphenol (catechol) moiety on a six-membered aromatic or conjugated ring"
    
    # If no matching ring was found, return a negative result.
    return False, "No adjacent qualifying hydroxyl substituents found on a six-membered aromatic or conjugated ring"

# Example usage (for testing):
if __name__ == '__main__':
    test_examples = [
        # True positives (catechols):
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
        ("C=1(C=CC=CC1C(=O)/C=C/C2=CC(=C(C=C2)O)O)O", "2',3,4-trihydroxy-trans-chalcone"),
        ("CN(C)CCC1=CC=CC=C1\\N=N\\C1=CC=C(O)C(O)=C1", "trans-azodopa"),
        ("CCCCCCCCCCCCCc1cc(C)c(C)c(O)c1O", "5,6-dimethyl-3-pentadecylcatechol"),
        # False positive examples (should not be classified as simple catechols)
        ("Oc1ccc(cc1OS([O-])(=O)=O)-c1oc2cc([O-])cc(O)c2c(=O)c1OS([O-])(=O)=O", "quercetin 3,3'-bissulfate(3-)"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_catechols(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {result} -> {reason}\n{'-'*80}")