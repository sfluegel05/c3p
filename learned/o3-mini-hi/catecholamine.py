"""
Classifies: CHEBI:33567 catecholamine
"""
#!/usr/bin/env python
"""
Classifies: catecholamine
Definition: Catecholamines include 4-(2-Aminoethyl)pyrocatechol and derivatives.
A valid catecholamine must contain an isolated aromatic benzene ring (six-membered),
with two adjacent oxygen substituents (which may be substituted variants),
and one ring atom must have a directly attached aminoethyl chain – that is, a chain
of exactly two aliphatic (sp³) carbon atoms connected by single bonds ending in a nitrogen.
This program uses RDKit for molecule parsing and heuristic substructure detection.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.

    Requirements:
      (i) An isolated aromatic benzene ring (six-membered) that has two adjacent oxygen substituents.
          The oxygen substituents (which may bear additional groups) must be connected by a single bond.
      (ii) One of the benzene ring atoms has a substituent chain of exactly length 3:
           the chain is: ring_atom -- (aliphatic/sp³ C) -- (aliphatic/sp³ C) -- (N)
           where the bonds along the chain are all single bonds.

    Args:
      smiles (str): A SMILES string representing the molecule.
    
    Returns:
      bool: True if the molecule is classified as a catecholamine, False otherwise.
      str: A textual reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get information about rings in the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Helper: check if an oxygen substituent is attached directly to an atom.
    def has_oxygen_substituent(atom, ring_atom_idxs):
        for nbr in atom.GetNeighbors():
            # only consider substituents outside the ring
            if nbr.GetIdx() in ring_atom_idxs:
                continue
            if nbr.GetAtomicNum() == 8:
                # allow oxygen if the bond is single (thus not part of a double bond carbonyl, etc)
                for b in mol.GetBonds():
                    # (not the most efficient but sufficient for our purposes)
                    pass
                # Check that the bond between atom and nbr is single.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    return True
        return False

    # Helper: check if two indices in a ring (cyclic order) are adjacent.
    # We compare each pair of indices that are consecutive in the sorted cyclic order.
    def adjacent_in_ring(ring, idx1, idx2):
        ring_list = list(ring)
        n = len(ring_list)
        for i in range(n):
            if (ring_list[i] == idx1 and ring_list[(i+1)%n] == idx2) or (ring_list[i] == idx2 and ring_list[(i+1)%n] == idx1):
                return True
        return False

    # Helper: check for an aminoethyl chain attached directly to a ring atom.
    # We expect a path: ring_atom -- C1 -- C2 -- N (each bond being single)
    # The two carbons must be aliphatic (non-aromatic) and their hybridization should be SP3.
    def has_aminoethyl_chain(ring_atom):
        for nbr1 in ring_atom.GetNeighbors():
            # only consider substituents not part of the ring 
            if nbr1.GetIsAromatic():
                continue
            # The bond from ring_atom to nbr1 should be single.
            b1 = mol.GetBondBetweenAtoms(ring_atom.GetIdx(), nbr1.GetIdx())
            if b1 is None or b1.GetBondType() != Chem.BondType.SINGLE:
                continue
            if nbr1.GetAtomicNum() != 6:
                continue
            # Check that nbr1 is sp³ hybridized
            if nbr1.GetHybridization() != rdchem.HybridizationType.SP3:
                continue
            # Look at neighbors of nbr1, but avoid going back to the ring_atom.
            for nbr2 in nbr1.GetNeighbors():
                if nbr2.GetIdx() == ring_atom.GetIdx():
                    continue
                b2 = mol.GetBondBetweenAtoms(nbr1.GetIdx(), nbr2.GetIdx())
                if b2 is None or b2.GetBondType() != Chem.BondType.SINGLE:
                    continue
                if nbr2.GetAtomicNum() != 6:
                    continue
                if nbr2.GetIsAromatic():
                    continue
                if nbr2.GetHybridization() != rdchem.HybridizationType.SP3:
                    continue
                # Now from nbr2, look for a nitrogen (any type) connected by a single bond.
                for nbr3 in nbr2.GetNeighbors():
                    if nbr3.GetIdx() == nbr1.GetIdx():
                        continue
                    b3 = mol.GetBondBetweenAtoms(nbr2.GetIdx(), nbr3.GetIdx())
                    if b3 is None or b3.GetBondType() != Chem.BondType.SINGLE:
                        continue
                    if nbr3.GetAtomicNum() == 7:
                        return True
        return False

    # Iterate over rings to find an isolated aromatic benzene ring.
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # get the atoms in this ring
        ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
        # All atoms in the ring must be aromatic.
        if not all(a.GetIsAromatic() for a in ring_atoms):
            continue
        # Exclude fused rings: each atom in the chosen ring should belong only once.
        is_isolated = True
        for idx in ring:
            count = sum(1 for r in ring_info if idx in r)
            if count > 1:
                is_isolated = False
                break
        if not is_isolated:
            continue
        
        ring_set = set(ring)
        # Find ring atom indices that have an oxygen substituent directly attached.
        oxy_idxs = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if has_oxygen_substituent(atom, ring_set):
                oxy_idxs.append(idx)
        if len(oxy_idxs) < 2:
            continue

        # Check whether there is at least one pair of adjacent ring atoms (neighbors within the ring)
        # that both have oxygen substituents.
        catechol_found = False
        for i in range(len(oxy_idxs)):
            for j in range(i+1, len(oxy_idxs)):
                if adjacent_in_ring(ring, oxy_idxs[i], oxy_idxs[j]):
                    catechol_found = True
                    break
            if catechol_found:
                break
        if not catechol_found:
            continue

        # Now check for a directly attached aminoethyl chain on any ring atom.
        if not any(has_aminoethyl_chain(mol.GetAtomWithIdx(idx)) for idx in ring):
            return False, "Catechol ring detected but no directly attached aminoethyl chain (-CH2-CH2-N) found"
        
        # If we reach here for at least one candidate ring then we classify this as a catecholamine.
        return True, "Contains an isolated catechol ring (with two adjacent oxygen substituents) and a directly attached aminoethyl chain"

    return False, "No isolated catechol ring (benzene with two adjacent oxygen substituents) found"


# For testing the function when run as a script.
if __name__ == "__main__":
    # Some example SMILES from the provided list.
    test_molecules = [
        ("C(CNCCCCCCNCCC1=CC=CC=C1)C2=CC(O)=C(C=C2)O", "dopexamine"),
        ("OC(=O)C1CC(=C\\C=N/CCc2ccc(O)c(O)c2)/C=C(N1)C(O)=O", "Miraxanthin-V"),
        ("C[C@H](N)[C@H](O)c1ccc(O)c(O)c1", "(-)-alpha-Methylnoradrenaline"),
        ("C=1(C=C(C(O)=CC1)O)CCN.Cl", "Dopamine hydrochloride"),
        ("C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O", "(R)-dobutamine"),
        ("C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1", "(S)-dobutamine"),
        ("OC(=O)[C@H](Cc1ccc(O)c(O)c1)\\N=C/C=C1C[C@H](NC(=C\\1)C(O)=O)C(O)=O", "Dopaxanthin"),
        ("CNC[C@H](O)c1ccc(O)c(O)c1", "(R)-adrenaline"),
        ("CC(N)C(O)c1ccc(O)c(O)c1", "Nordephrine"),
        ("NC[C@@H](O)c1ccc(O)c(O)c1", "(S)-noradrenaline"),
        ("[C@@H]([C@@H](N)C)(O)C1=CC(O)=C(C=C1)O", "alpha-methylnoradrenaline"),
        ("C1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN", "4-(2-aminoethyl)-5-nitrobenzene-1,2-diol"),
        ("C1=CC(=C(C(=C1CCN)O)O)O", "4-(2-aminoethyl)benzene-1,2,3-triol"),
        ("O[C@@H](CNCCCCc1ccc(O)cc1)c1ccc(O)c(O)c1", "arbutamine"),
        ("CNC[C@@H](O)c1ccc(O)c(O)c1", "(S)-adrenaline"),
        ("CNCC(O)C1=CC(O)=C(O)C=C1", "4-[1-hydroxy-2-(methylamino)ethyl]benzene-1,2-diol"),
        ("NC[C@H](O)c1ccc(O)c(O)c1", "(R)-noradrenaline"),
        ("CC(C)NC[C@H](O)c1ccc(O)c(O)c1", "L-isoprenaline"),
        ("NCC(O)c1ccc(O)c(O)c1", "noradrenaline"),
        ("OC1=C(O)C=C(CCN)C=C1O", "5-Hydroxydopamine"),
        ("N[C@@H](CSc1cc(C[C@H](N)C(O)=O)cc(O)c1O)C(O)=O", "Cysteinyldopa"),
        ("CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1", "dobutamine"),
        ("NCCc1ccc(O)c(O)c1", "dopamine"),
        ("CCCN(CCC)CCC1=CC(=C(C=C1)O)O", "4-[2-(dipropylamino)ethyl]benzene-1,2-diol"),
        ("C[N+](C)(C)CCc1ccc(O)c(O)c1", "Coryneine")
    ]
    
    for smi, name in test_molecules:
        result, reason = is_catecholamine(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")