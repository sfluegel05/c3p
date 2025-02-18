"""
Classifies: CHEBI:33567 catecholamine
"""
#!/usr/bin/env python
"""
Classifies: catecholamine
Definition: Catecholamines include 4-(2-Aminoethyl)pyrocatechol and derivatives.
A valid catecholamine must contain:
  (i) a catechol moiety – that is an aromatic six-membered ring (benzene) with two adjacent oxygen substituents (e.g. hydroxyl groups),
  (ii) an aminoethyl side chain – a 2-carbon aliphatic chain ending in a nitrogen – that is attached directly to one of the ring atoms.
This program uses RDKit for molecule parsing and heuristic substructure detection.
"""

from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    It requires:
      - An aromatic six-membered ring that has two adjacent oxygen substituents (hydroxyl-like group) attached,
      - An aminoethyl side chain (CH2-CH2-N...) attached directly to an atom of that catechol ring.
    
    Args:
      smiles (str): SMILES string representing the molecule.
    
    Returns:
      bool: True if the molecule is classified as a catecholamine, False otherwise.
      str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo().AtomRings()
    found_catechol_ring = False
    found_aminoethyl = False
    catechol_ring_indices = None

    # Helper function: Check if a given oxygen atom attached to a ring atom looks like a hydroxyl.
    def is_hydroxyl(o_atom):
        # Check oxygen: must be atomic number 8, attached via a single bond,
        # and have at least one hydrogen (explicit or implicit).
        if o_atom.GetAtomicNum() != 8:
            return False
        # Must be connected via a SINGLE bond.
        for bond in o_atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False
        # Check hydrogen count: note that RDKit usually has implicit Hs.
        if o_atom.GetTotalNumHs() < 1:
            return False
        return True

    # Helper function: Check if a substituent on a ring atom is an aminoethyl chain.
    def has_aminoethyl_chain(ring_atom, ring_set):
        # Look at each neighbor of ring_atom that is not in the ring.
        for nbr1 in ring_atom.GetNeighbors():
            if nbr1.GetIdx() in ring_set:
                continue
            # The first atom should be a non‐aromatic carbon (ideally CH2, but may have extra substituents)
            if nbr1.GetAtomicNum() != 6 or nbr1.GetIsAromatic():
                continue
            # Now, look for a second carbon attached to nbr1 (and not part of the ring and not going back to ring_atom)
            for nbr2 in nbr1.GetNeighbors():
                if nbr2.GetIdx() == ring_atom.GetIdx():
                    continue
                if nbr2.GetIdx() in ring_set:
                    continue
                if nbr2.GetAtomicNum() != 6 or nbr2.GetIsAromatic():
                    continue
                # Finally, check if nbr2 has a nitrogen attached (other than nbr1)
                for nbr3 in nbr2.GetNeighbors():
                    if nbr3.GetIdx() in (nbr1.GetIdx(), ring_atom.GetIdx()):
                        continue
                    if nbr3.GetAtomicNum() == 7:
                        return True
        return False

    # Loop over rings – we only care about aromatic six-membered rings.
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Check that every atom of the ring is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # Now check for catechol: find two adjacent ring atoms that both have an oxygen substituent (likely –OH).
        # Instead of using cyclic order from the ring definition (which may be in arbitrary order),
        # we check each bond within the ring.
        ring_set = set(ring)
        # For each atom in the ring, note if it has a hydroxyl substituent attached.
        hydroxyl_on_atom = {}
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            has_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                if is_hydroxyl(nbr):
                    has_OH = True
                    break
            hydroxyl_on_atom[idx] = has_OH

        # Now check bonds between atoms in the ring: if two bonded ring atoms both have OH.
        catechol_in_ring = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                if a1.GetIdx() in ring_set and a2.GetIdx() in ring_set:
                    if hydroxyl_on_atom.get(a1.GetIdx(), False) and hydroxyl_on_atom.get(a2.GetIdx(), False):
                        catechol_in_ring = True
                        break
            if catechol_in_ring:
                break

        if not catechol_in_ring:
            continue

        # We have a catechol ring. Now check if at least one ring atom has an attached aminoethyl chain.
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            if has_aminoethyl_chain(ring_atom, ring_set):
                found_aminoethyl = True
                break

        if catechol_in_ring:
            found_catechol_ring = True
            catechol_ring_indices = ring_set
            # If we found an aminoethyl chain for this ring, we can stop.
            if found_aminoethyl:
                break

    if not found_catechol_ring:
        return False, "No catechol ring (aromatic six‐membered ring with two adjacent hydroxyl groups) found"

    if not found_aminoethyl:
        return False, "Catechol ring detected but no directly attached aminoethyl chain (–CH2–CH2–N...) found"

    return True, "Contains a catechol ring with an aminoethyl side chain attached directly to the ring"

# For testing purposes when run as a stand-alone script.
if __name__ == "__main__":
    test_molecules = [
        # True positives (catecholamine examples)
        ("C(CNCCCCCCNCCC1=CC=CC=C1)C2=CC(O)=C(C=C2)O", "dopexamine"),
        ("OC(=O)C1CC(=C\\C=N/CCc2ccc(O)c(O)c2)/C=C(N1)C(O)=O", "Miraxanthin-V"),
        ("C=1(C=C(C(=CC1)O)OC)C(O)CN", "Normetanephrine"),
        ("C[C@H](N)[C@H](O)c1ccc(O)c(O)c1", "(-)-alpha-Methylnoradrenaline"),
        ("C=1(C=C(C(O)=CC1)O)CCN.Cl", "Dopamine hydrochloride"),
        ("C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O", "(R)-dobutamine"),
        ("C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1", "(S)-dobutamine"),
        ("OC(=O)[C@H](Cc1ccc(O)c(O)c1)\\N=C/C=C1C[C@H](NC(=C\\1)C(O)=O)C(O)=O", "Dopaxanthin"),
        ("CNC[C@H](O)c1ccc(O)c(O)c1", "(R)-adrenaline"),
        ("CC(N)C(O)c1ccc(O)c(O)c1", "Nordephrine"),
        ("O(C)C1=C(O)C=CC(C(O)CNC)=C1", "Metanephrine"),
        ("NC[C@@H](O)c1ccc(O)c(O)c1", "(S)-noradrenaline"),
        ("[C@@H]([C@@H](N)C)(O)C1=CC(O)=C(C=C1)O", "alpha-methylnoradrenaline"),
        ("C1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN", "4-(2-aminoethyl)-5-nitrobenzene-1,2-diol"),
        ("C1=CC(=C(C(=C1CCN)O)O)O", "4-(2-aminoethyl)benzene-1,2,3-triol"),
        ("O[C@@H](CNCCCCc1ccc(O)cc1)c1ccc(O)c(O)c1", "arbutamine"),
        ("CNC[C@@H](O)c1ccc(O)c(O)c1", "(S)-adrenaline"),
        ("CNCC(O)C1=CC(O)=C(O)C=C1", "4-[1-hydroxy-2-(methylamino)ethyl]benzene-1,2-diol"),
        ("NC[C@H](O)c1ccc(O)c(O)c1", "(R)-noradrenaline"),
        ("CC(C)NC[C@H](O)c1ccc(O)c(O)c1", "L-isoprenaline"),
        ("S(OC1=C(O)C=C([C@@H](O)CNC)C=C1)(O)(=O)=O", "Epinephrine sulfate"),
        ("NCC(O)c1ccc(O)c(O)c1", "noradrenaline"),
        ("OC1=C(O)C=C(CCN)C=C1O", "5-Hydroxydopamine"),
        ("N[C@@H](CSc1cc(C[C@H](N)C(O)=O)cc(O)c1O)C(O)=O", "Cysteinyldopa"),
        ("CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1", "dobutamine"),
        ("NCCc1ccc(O)c(O)c1", "dopamine"),
        ("CCCN(CCC)CCC1=CC(=C(C=C1)O)O", "4-[2-(dipropylamino)ethyl]benzene-1,2-diol"),
        ("C[N+](C)(C)CCc1ccc(O)c(O)c1", "Coryneine"),
        # False positives (molecules that have similar fragments but should not be classified as catecholamines)
        ("CCOC1=C(C=C(C=C1)CCNC(=O)CC2=CC=CC=C2)OCC", "N-[2-(3,4-diethoxyphenyl)ethyl]-2-phenylacetamide"),
        # Additional false positive examples provided can be tested similarly.
    ]
    
    for smi, name in test_molecules:
        result, reason = is_catecholamine(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")