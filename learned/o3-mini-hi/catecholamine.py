"""
Classifies: CHEBI:33567 catecholamine
"""
#!/usr/bin/env python
"""
Classifies: catecholamine
Definition: Catecholamines include 4-(2-Aminoethyl)pyrocatechol and derivatives.
A valid catecholamine must contain an isolated catechol ring – an aromatic benzene ring with two adjacent oxygen‐substituents (which may be –OH or substituted variants),
and one of the ring atoms must have a directly attached aminoethyl (–CH2–CH2–N…) chain.
This program uses RDKit for molecule parsing and heuristic substructure detection.
"""

from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    Requirements:
      (i) An isolated aromatic six‐membered ring (benzene) that has two adjacent oxygen substituents. 
          The oxygen need not be a free hydroxyl; substitutions (like –OCH3, –OSO3H) are allowed.
      (ii) A directly attached aminoethyl side chain, meaning that one ring atom has a single-bonded chain 
           that is exactly two non‐aromatic carbons long ending in a nitrogen atom.
    
    Args:
      smiles (str): SMILES string representing the molecule.
    
    Returns:
      bool: True if molecule is classified as a catecholamine, False otherwise.
      str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Helper: Determine if an oxygen substituent is attached directly.
    def is_oxygen_substituent(o_atom):
        # Accept oxygen if it is singly bonded (so it is not a carbonyl or in a delocalized system).
        # We do not insist on a hydrogen because it might be substituted.
        if o_atom.GetAtomicNum() != 8:
            return False
        for bond in o_atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False
        return True

    # Helper: Check if a given ring atom has an oxygen substituent (attached from outside the ring).
    def has_oxygen_substituent(atom, ring_set):
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in ring_set:
                continue
            if is_oxygen_substituent(nbr):
                return True
        return False

    # Helper: Check if the substituent on this ring atom forms an aminoethyl chain.
    # We require a chain: ring_atom - C1 - C2 - N (where C1 and C2 are non-aromatic sp3 carbons not belonging to any aromatic ring).
    def has_aminoethyl_chain(ring_atom):
        for nbr1 in ring_atom.GetNeighbors():
            # The neighbor should not be part of any aromatic ring and be carbon.
            if nbr1.GetAtomicNum() != 6 or nbr1.GetIsAromatic():
                continue
            # To ensure we have a chain (and not some cyclic structure), we can require nbr1 is not in any ring.
            if nbr1.IsInRing():
                continue
            # Now examine neighbors of nbr1 (except ring_atom)
            for nbr2 in nbr1.GetNeighbors():
                if nbr2.GetIdx() == ring_atom.GetIdx():
                    continue
                if nbr2.GetAtomicNum() != 6 or nbr2.GetIsAromatic():
                    continue
                if nbr2.IsInRing():
                    continue
                # Now from nbr2, look for at least one nitrogen neighbor (not going back to nbr1)
                for nbr3 in nbr2.GetNeighbors():
                    if nbr3.GetIdx() in (nbr1.GetIdx(), ring_atom.GetIdx()):
                        continue
                    if nbr3.GetAtomicNum() == 7:
                        return True
        return False

    # Now, iterate over candidate rings: isolated, aromatic benzene rings of size 6.
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is aromatic.
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if not all(atom.GetIsAromatic() for atom in atoms):
            continue

        # Exclude fused rings: require that each atom in the candidate ring occurs in only one ring.
        is_isolated = True
        for idx in ring:
            count = sum(1 for r in ring_info if idx in r)
            if count > 1:
                is_isolated = False
                break
        if not is_isolated:
            continue

        ring_set = set(ring)
        # Find all ring atoms that have an oxygen substituent.
        oxy_atoms = [idx for idx in ring if has_oxygen_substituent(mol.GetAtomWithIdx(idx), ring_set)]
        if len(oxy_atoms) < 2:
            continue

        # Now check that there exists at least one pair of adjacent atoms in the ring (adjacent via a ring bond)
        # that both have an oxygen substituent.
        catechol_found = False
        for idx in oxy_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set and nbr.GetIdx() in oxy_atoms:
                    catechol_found = True
                    break
            if catechol_found:
                break
        if not catechol_found:
            continue

        # At this point, we have a candidate catechol ring.
        # Now look for an attached aminoethyl chain from any ring atom.
        aminoethyl_found = False
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            if has_aminoethyl_chain(ring_atom):
                aminoethyl_found = True
                break
        
        if not aminoethyl_found:
            return False, "Catechol ring detected but no directly attached aminoethyl chain (-CH2-CH2-N...) found"
        
        # If we reach here for a candidate ring, we consider the molecule a catecholamine.
        return True, "Contains an isolated catechol ring (with two adjacent oxygen substituents) and a directly attached aminoethyl chain"
    
    return False, "No isolated catechol ring (benzene with two adjacent oxygen substituents) found"

# For testing purposes when run as a stand-alone script.
if __name__ == "__main__":
    test_molecules = [
        # True positives (catecholamine examples)
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
        ("C[N+](C)(C)CCc1ccc(O)c(O)c1", "Coryneine"),
        # False positives (molecules that should not be classified as catecholamines)
        ("Oc1ccc(C[C@@H]2[NH2+]CCc3cc(O)c(O)cc23)cc1", "(S)-norcoclaurinium(1+)"),
        ("O=C(O)[C@H]1N=CC=2C=C(O)C(=CC2C1)O", "JS-1"),
        ("OC1CN(CC=2C1=CC(O)=C(O)C2)C", "N-methyl-4,6,7-trihydroxy-1,2,3,4-tetrahydroisoquinoline"),
        ("OC1=C(O)C=C2[NH2+][C@@H](CC2=C1)C([O-])=O", "leucodopachrome zwitterion"),
        ("O=C(NCCC1=CC(O)=C(O)C=C1)CCC2=CC(O)=C(O)C=C2", "3-(3,4-dihydroxyphenyl)-N-[2-(3,4-dihydroxyphenyl)ethyl]propanimidic acid"),
        ("CC1=NC=CN1CC(=O)C2=CC(=C(C=C2)O)O", "1-(3,4-dihydroxyphenyl)-2-(2-ethyl-1-imidazolyl)ethanone"),
        ("Br.C1CNC[C@@H](C2=C1C=C(C(=C2)O)O)C3=CC=CC=C3", "(R)-SKF 38393 hydrobromide"),
    ]
    
    for smi, name in test_molecules:
        result, reason = is_catecholamine(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")