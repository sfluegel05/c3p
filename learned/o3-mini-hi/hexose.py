"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: hexose
Definition: Any six‐carbon monosaccharide which in its linear form contains either an aldehyde group
at position 1 (aldohexose) or a ketone group at position 2 (ketohexose). In its cyclic form the sugar
often forms a pyranose (6–membered ring, 5 C + 1 O plus an exocyclic –CH2OH) or a furanose (5–membered ring
plus an exocyclic CH2OH) ring. (All criteria are heuristic.)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule qualifies as a hexose.
    
    The strategy is three–fold:
      1. Quick rejection: if the molecule contains any phosphorus.
      2. Open–chain (linear) detection: using relaxed SMARTS for aldo–hexose and keto–hexose.
         (We do not insist on explicit H counts so as to allow small substitutions.)
      3. Cyclic (ring) detection: look through all rings for a pyranose candidate (ring of 6 atoms with exactly one O and 5 C)
         or a furanose candidate (ring of 5 atoms with exactly one O plus exocyclic carbon(s) that bring the total sugar C count to 6).
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a hexose, otherwise False.
        str: Explanation for the classification decision.
    """
    # Parse SMILES and add hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # 1. Quick rejection: if molecule contains phosphorus (often a nucleotide or other non‐sugar)
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Contains phosphorus, likely not a hexose"
    
    # ------------------------
    # 2. Open–chain detection using relaxed SMARTS.
    # For an aldo–hexose, the typical pattern is: O=CH-(CHOH)4-CH2OH.
    # We use a relaxed SMARTS (dropping explicit H counts) so that minor substitutions still match.
    aldo_smarts = ("[CX3](=O)"                # C1: aldehyde group
                   "[CX4]([OX2])"             # C2: CH–OH division (allow substitution on –OH)
                   "[CX4]([OX2])"             # C3
                   "[CX4]([OX2])"             # C4
                   "[CX4]([OX2])"             # C5
                   "[CX4][OX2]")              # C6: CH2–OH; note we do not force [CX4H2] here
    aldo_query = Chem.MolFromSmarts(aldo_smarts)
    if aldo_query is None:
        return None, None  # Should not happen
    matches = mol.GetSubstructMatches(aldo_query)
    for match in matches:
        # Count how many carbons are in the matching fragment.
        n_carbons = sum(1 for idx in match if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
        if n_carbons == 6:
            return True, "Matches open–chain aldo–hexose pattern (aldehyde at C1)"
    
    # For a keto–hexose (typical open–chain structure: CH2OH–C(=O)–(CHOH)3–CH2OH)
    keto_smarts = ("[CX4][OX2]"              # C1: CH2OH
                   "[CX3](=O)"               # C2: ketone group
                   "[CX4]([OX2])"            # C3: CH–OH
                   "[CX4]([OX2])"            # C4: CH–OH
                   "[CX4]([OX2])"            # C5: CH–OH
                   "[CX4][OX2]")             # C6: CH2OH
    keto_query = Chem.MolFromSmarts(keto_smarts)
    if keto_query is None:
        return None, None
    matches = mol.GetSubstructMatches(keto_query)
    for match in matches:
        n_carbons = sum(1 for idx in match if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
        if n_carbons == 6:
            return True, "Matches open–chain keto–hexose pattern (ketone at C2)"
    
    # ------------------------
    # 3. Cyclic (ring) detection
    # The idea is to search for rings that seem to be sugar cores.
    # For pyranose: look for a ring of 6 atoms with 1 oxygen and 5 carbons. Then look for at least one exocyclic carbon
    # attached to the ring carbons (most often CH2OH) to bring the total sugar C count to 6.
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    for ring in rings:
        if len(ring) == 6:
            # pyranose candidate: expect 5 carbons and 1 oxygen in the ring
            ring_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            n_carbon = ring_symbols.count("C")
            n_oxygen = ring_symbols.count("O")
            if n_carbon == 5 and n_oxygen == 1:
                # Look for an exocyclic carbon attached to any ring carbon.
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() != "C":
                        continue
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in ring:
                            continue
                        if nbr.GetSymbol() == "C":
                            # For our purposes, check that this substituent carbon has at least one oxygen neighbor.
                            if any(n.GetSymbol() == "O" for n in nbr.GetNeighbors() if n.GetIdx() != idx):
                                # We now have a pyranose core of 5 ring carbons plus one exocyclic carbon.
                                return True, "Contains pyranose ring pattern (6-membered ring with one oxygen and exocyclic CH2OH branch)"
        # For a furanose candidate: look for a ring of 5 atoms.
        if len(ring) == 5:
            # furanose candidate: expect 4 carbons and 1 oxygen in ring.
            ring_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            n_carbon = ring_symbols.count("C")
            n_oxygen = ring_symbols.count("O")
            if n_carbon == 4 and n_oxygen == 1:
                # For many hexoses in furanose form, the sugar core is made up of the 4 ring carbons plus 2 exocyclic carbons.
                exo_carbons = set()
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetSymbol() != "C":
                        continue
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() in ring:
                            continue
                        if nbr.GetSymbol() == "C":
                            exo_carbons.add(nbr.GetIdx())
                if (n_carbon + len(exo_carbons)) == 6:
                    return True, "Contains furanose ring pattern (5-membered ring with one oxygen and sufficient exocyclic carbons)"
    
    # If none of the above checks are satisfied, then we do not classify the molecule as hexose.
    return False, "Does not match recognized hexose patterns"

# ----- Example usage -----
if __name__ == '__main__':
    # Here are some test cases (both positive and problematic cases):
    test_smiles = [
        "OC(C(O)CNCCCCCCC)C(O)C(O)CO",                         # 1-Deoxy-1-(heptylamino)hexitol (modified open-chain hexose)
        "O=C(OC1OC(C(O)C(C1O)O)C)C2=CC=CC=C2",                 # 1-O-Benzoyl-alpha-L-rhamnopyranoside (cyclic, glycoside)
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",               # D-allopyranose (cyclic pyranose)
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",         # aldehydo-D-galactose (open-chain aldo)
        "C[C@H](O)[C@H](O)[C@@H](O)C(=O)CO",                   # L-rhamnulose (open-chain keto variant)
        "[H]C([H])([C@]([H])(O)C=O)[C@]([H])(O)[C@@]([H])(C)O",  # tyvelose (open-chain with deoxy substitution)
        "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H]1O",               # alpha-D-tagatofuranose (cyclic furanose)
    ]
    for s in test_smiles:
        result, reason = is_hexose(s)
        print(f"SMILES: {s}\n  Classified as hexose? {result}. Reason: {reason}\n")