"""
Classifies: CHEBI:33567 catecholamine
"""
#!/usr/bin/env python
"""
Classifies: catecholamine (4-(2-Aminoethyl)pyrocatechol and derivatives)
A catecholamine is defined by the presence of a catechol moiety (an aromatic benzene ring
bearing two adjacent oxygen substituents, which may be hydroxyl or substituted forms)
and an aminoethyl side chain (a two‐carbon chain terminating in a nitrogen, allowing common substitutions).
This program uses heuristic substructure detection to classify a SMILES string.
"""

from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    The heuristic rules are:
      1. The molecule must contain at least one aromatic benzene ring (6-membered) 
         that features two adjacent ring atoms each substituted (via a single bond) with an oxygen.
         (This oxygen need not be a free –OH but may be substituted (e.g. methoxy or sulfate esters).)
      2. The molecule must contain an aminoethyl side chain attached directly to the aromatic ring.
         Since variations occur (e.g. a hydroxyl substitution at the carbon alpha to the ring),
         several SMARTS patterns are examined that cover common variations.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a catecholamine, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- STEP 1: Detect a catechol (benzene-derived) moiety ---
    # We look over 6-membered rings that are aromatic. Then we check if at least two adjacent
    # ring atoms carry an oxygen substituent via a SINGLE bond.
    catechol_found = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue
        # Ensure all atoms of the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue

        # Identify which ring atoms have a single-bonded oxygen neighbor (allowing for -OH, -OCH3, -OSO3, etc).
        oxy_substituted = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Exclude atoms that are in the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Check if the neighbor is oxygen and the bond is a single bond.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if nbr.GetAtomicNum() == 8 and bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    oxy_substituted.add(idx)
                    break  # Once found one oxygen, move to next ring atom.
        # Now check if there are at least two adjacent atoms (in cyclic order) among the substituted ones.
        if len(oxy_substituted) < 2:
            continue
        ring_len = len(ring)
        for i in range(ring_len):
            current_atom = ring[i]
            next_atom = ring[(i+1) % ring_len]
            if current_atom in oxy_substituted and next_atom in oxy_substituted:
                catechol_found = True
                break
        if catechol_found:
            break
    
    if not catechol_found:
        return False, "No catechol moiety (aromatic ring with two adjacent oxygen substituents) found"
    
    # --- STEP 2: Detect an aminoethyl side chain ---
    # We expect a two-carbon chain from an aromatic carbon leading to an amine.
    # However, in derivatives the first carbon may be substituted (e.g. carrying an -OH).
    # We test several SMARTS patterns to cover common variants.
    aminoethyl_smarts = [
        "cCCN",         # simple: aromatic carbon - CH2 - CH2 - N(any substitution)
        "cC(O)CN",      # variant: aromatic carbon - CH(OH) - CH2 - N
        "cC([OX2H])CN", # alternative explicit pattern for free hydroxyl on the first carbon.
        "c[C@H](O)CN",  # chiral variant
        "c[C@@H](O)CN"  # chiral variant
    ]
    
    side_chain_found = False
    for smarts in aminoethyl_smarts:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            side_chain_found = True
            break
    
    if not side_chain_found:
        return False, "No aminoethyl side chain (e.g. pattern resembling 'cCCN' or derivatives) found"
    
    # If both key features are present, we classify as catecholamine.
    return True, "Contains a catechol moiety with an aminoethyl side chain"

# Example test cases
if __name__ == "__main__":
    test_data = [
        # True positives
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
        # False negatives (expected to be matched, so if not, we output details):
        ("C=1(C=C(C(=CC1)O)OC)C(O)CN", "Normetanephrine"),
        ("S(OC1=C(O)C=C([C@@H](O)CN)C=C1)(O)(=O)=O", "Norepinephrine sulfate"),
        ("O(C)C1=C(O)C=CC(C(O)CNC)=C1", "Metanephrine"),
        ("S(OC1=C(O)C=C([C@@H](O)CNC)C=C1)(O)(=O)=O", "Epinephrine sulfate")
    ]
    
    for smi, name in test_data:
        result, reason = is_catecholamine(smi)
        print(f"SMILES: {smi}\nName: {name}\nResult: {result}\nReason: {reason}\n{'-'*60}")