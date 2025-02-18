"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: Monoamine
Definition: A monoamine is an aralkylamino compound having exactly one amino group that is connected 
to an aromatic ring by a two‐carbon chain. In this implementation we “walk” the molecular graph: 
for each non‐aromatic, non‐ring nitrogen (that is not directly attached to a carbonyl) we search for a 
unique path: N – C1 – C2 – aromatic C (all bonds single; C1 and C2 are aliphatic).
If exactly one such aminic moiety exists then the molecule is classified as a monoamine.
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.

    A monoamine is defined as an aralkylamino compound having exactly one amino group connected to an 
    aromatic ring by a two‐carbon chain. To decide this we traverse the molecular graph. For each nitrogen
    atom that is (a) non-aromatic, (b) not in a ring, and (c) not directly attached to a carbonyl carbon,
    we try to find a single path of three bonds: from that nitrogen to an aliphatic carbon (C1), then to a 
    second aliphatic carbon (C2) and finally to an aromatic carbon. Only single bonds are allowed 
    along the chain and the intermediary carbons (C1 and C2) must be non-aromatic and not in a ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains exactly one monoamine moiety, False otherwise.
        str: A message explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Helper function: returns True if the given atom is a carbonyl carbon (i.e. bonded
    # by a double bond to an oxygen) – used to disqualify N atoms directly attached to a carbonyl.
    def is_adjacent_to_carbonyl(atom):
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon
                for bond in nbr.GetBonds():
                    # Check other atom of the bond
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        return True
        return False

    # This function tries to walk from candidate N to an aromatic carbon through exactly 2 carbons.
    def has_valid_two_carbon_link(n_atom):
        # Examine neighbors of n_atom that are carbon and satisfy conditions:
        for c1 in n_atom.GetNeighbors():
            if c1.GetAtomicNum() != 6: 
                continue
            # c1 must be non-aromatic, not in a ring and bond must be SINGLE.
            bond1 = mol.GetBondBetweenAtoms(n_atom.GetIdx(), c1.GetIdx())
            if bond1.GetBondType() != Chem.BondType.SINGLE:
                continue
            if c1.GetIsAromatic() or c1.IsInRing():
                continue
            # Now from c1, find a neighboring carbon (C2) different from n_atom:
            for c2 in c1.GetNeighbors():
                if c2.GetIdx() == n_atom.GetIdx():
                    continue
                if c2.GetAtomicNum() != 6:
                    continue
                bond2 = mol.GetBondBetweenAtoms(c1.GetIdx(), c2.GetIdx())
                if bond2.GetBondType() != Chem.BondType.SINGLE:
                    continue
                if c2.GetIsAromatic() or c2.IsInRing():
                    continue
                # Finally, from c2, look for an aromatic carbon (A) not equal to c1:
                for a in c2.GetNeighbors():
                    if a.GetIdx() == c1.GetIdx():
                        continue
                    if a.GetAtomicNum() == 6 and a.GetIsAromatic():
                        # Found a valid N -> C1 -> C2 -> aromatic C path.
                        return True
        return False

    monoamine_nitrogens = set()
    # Iterate over all atoms and select candidate nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # Only consider non-aromatic nitrogen not residing in a ring.
        if atom.GetIsAromatic() or atom.IsInRing():
            continue
        # Exclude nitrogen atoms attached to a carbonyl carbon.
        if is_adjacent_to_carbonyl(atom):
            continue
        # Now, does this nitrogen have a valid path to an aromatic ring via a 2-carbon chain?
        if has_valid_two_carbon_link(atom):
            monoamine_nitrogens.add(atom.GetIdx())
    
    count = len(monoamine_nitrogens)

    if count == 0:
        return False, "No aralkylamino moiety (amine connected by a two‐carbon chain to an aromatic ring) found"
    elif count > 1:
        return False, f"Found {count} aralkylamino moieties; expected exactly one monoamine group"
    else:
        return True, "Contains exactly one amino group connected via a two‐carbon chain to an aromatic ring"

# Example usage (for manual testing):
if __name__ == "__main__":
    # Test examples (from given outcomes):
    test_smiles = [
        "CNC[C@@H](O)c1ccc(O)c(O)c1", # (S)-adrenaline: expected True.
        "NCCc1ccc(O)cc1",             # tyramine: expected True.
        "OC(=O)[C@H](Cc1ccc(O)c(O)c1)\\N=C/C=C1C[C@H](NC(=C\\1)C(O)=O)C(O)=O", # Dopaxanthin: expected True.
        "N[C@@H](CSc1cc(C[C@H](N)C(O)=O)cc(O)c1O)C(O)=O", # Cysteinyldopa: True.
        "CC(C)NC[C@H](O)c1ccc(O)c(O)c1",  # L-isoprenaline: True.
        "NC[C@H](O)c1ccc(O)c(O)c1",       # (R)-noradrenaline: True.
        "C=1(C=C(C(O)=CC1)O)CCN.Cl",       # Dopamine hydrochloride: True.
        "OC(=O)C1CC(=C\\C=N/CCc2ccc(O)c(O)c2)/C=C(N1)C(O)=O", # Miraxanthin-V: True.
        "C=1(C=C(C(=CC1)O)OC)C(O)CN",      # Normetanephrine: True.
        "O[C@@H](CNCCCCc1ccc(O)cc1)c1ccc(O)c(O)c1", # arbutamine: True.
        "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",    # (S)-dobutamine: True.
        "CNC[C@H](O)c1ccc(O)c(O)c1",      # (R)-adrenaline: True.
        "C=1(C(=CC=C(C1)CCN[C@@H](CCC=2C=CC(=CC2)O)C)O)O", # (R)-dobutamine: True.
        "O(C)C1=C(O)C=CC(C(O)CNC)=C1",    # Metanephrine: True.
        "CNCC(O)C1=CC(O)=C(O)C=C1",       # 4-[1-hydroxy-2-(methylamino)ethyl]benzene-1,2-diol: True.
        "C[N+](C)(C)CCc1ccc(O)c(O)c1",    # Coryneine: True.
        "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",# dobutamine: True.
        "[C@@H]([C@@H](N)C)(O)C1=CC(O)=C(C=C1)O", # alpha-methylnoradrenaline: True.
        "[NH3+]CCc1ccc(O)cc1",            # tyraminium: True.
        "C1=CC(=C(C(=C1CCN)O)O)O",         # 4-(2-aminoethyl)benzene-1,2,3-triol: True.
        "NCC(O)c1ccc(O)c(O)c1",           # noradrenaline: True.
        "CCCN(CCC)CCC1=CC(=C(C=C1)O)O",    # 4-[2-(dipropylamino)ethyl]benzene-1,2-diol: True.
        "C1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN", # 4-(2-aminoethyl)-5-nitrobenzene-1,2-diol: True.
        "OC1=C(O)C=C(CCN)C=C1O",          # 5-Hydroxydopamine: True.
        "C[C@H](N)[C@H](O)c1ccc(O)c(O)c1", # (-)-alpha-Methylnoradrenaline: True.
        "NCCc1ccc(O)c(O)c1",             # dopamine: True.
        "S(OC1=C(O)C=C([C@@H](O)CNC)C=C1)(O)(=O)=O", # Epinephrine sulfate: True.
        "CC(N)C(O)c1ccc(O)c(O)c1",        # Nordephrine: True.
        "NC[C@@H](O)c1ccc(O)c(O)c1",       # (S)-noradrenaline (again): True.
        "CC1=C(C(=CC(=C1O)O)O)CCN",        # 5-(2-aminoethyl)-6-methylbenzene-1,2,4-triol: True.
        "S(OC1=C(O)C=C([C@@H](O)CN)C=C1)(O)(=O)=O"  # Norepinephrine sulfate: True.
    ]
    for s in test_smiles:
        result, reason = is_monoamine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")