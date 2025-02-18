"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: hexose
Definition: Any six‐carbon monosaccharide which in its linear form contains either an aldehyde group
at position 1 (aldohexose) or a ketone group at position 2 (ketohexose). In its cyclic form the sugar
often forms a pyranose (6‐membered ring, 5 C + 1 O plus an exocyclic –CH2OH) or a furanose (5–membered ring
plus an exocyclic CH2OH) ring. The algorithm below uses a combination of strict SMARTS for open–chain forms
and ring‐analysis for cyclic forms. (All criteria are heuristic.)
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule qualifies as a hexose.
    
    The strategy is three–fold:
      1. Quick rejection: if the molecule contains phosphorus.
      2. Open–chain detection: using tight SMARTS for an aldo–hexose and a keto–hexose.
         These SMARTS explicitly require a six–carbon backbone decorated with the expected hydroxyl groups.
      3. Cyclic (ring) detection: look through all rings. For a pyranose candidate, we expect a 6–membered ring
         with exactly one oxygen (i.e. 5 C in the ring) plus a single exocyclic CH2OH branch on one of the ring carbons.
         For a furanose candidate, we expect a 5–membered ring (4 C, 1 O) plus an exocyclic CH2OH branch.
         In either case we also require that the candidate substructure is not just a minor decoration (its heavy–atom count
         must be at least about 30% of the whole molecule).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a hexose, otherwise False.
        str: Explanation for the classification decision.
    """
    # Parse SMILES and add implicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # (Working with explicit Hs may help in counting –OH groups)
    mol = Chem.AddHs(mol)
    
    # 1. Quick rejection: If molecule contains phosphorus, assume it is not a hexose.
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Contains phosphorus, likely a nucleotide or similar structure"
    
    # -------------------------
    # 2. Open–chain detection using SMARTS.
    #
    # For aldo–hexose (chain numbering: C1=CHO, C2–C5: each with –OH, C6: CH2OH)
    # The SMARTS below is built so that branching (the –OH groups) are explicitly matched.
    # Each bracketed group (e.g. "[CX3H1](=O)") counts as one atom, and the branch in parentheses
    # means the oxygen is matched as a separate atom.
    aldo_smarts = ("[CX3H1](=O)"               # C1: aldehyde (CHO)
                   "[CX4H]([OX2H])"            # C2: CH–OH
                   "[CX4H]([OX2H])"            # C3: CH–OH
                   "[CX4H]([OX2H])"            # C4: CH–OH
                   "[CX4H]([OX2H])"            # C5: CH–OH
                   "[CX4H2][OX2H]")            # C6: CH2–OH
    aldo_query = Chem.MolFromSmarts(aldo_smarts)
    if aldo_query is None:
        # Should not happen, but just in case.
        return None, None
    # For keto–hexose (typical open–chain structure: CH2OH–C(=O)–CHOH–CHOH–CHOH–CH2OH)
    keto_smarts = ("[CX4H2][OX2H]"           # C1: CH2–OH
                   "[CX3](=O)"                # C2: C=O (ketone)
                   "[CX4H]([OX2H])"           # C3: CH–OH
                   "[CX4H]([OX2H])"           # C4: CH–OH
                   "[CX4H]([OX2H])"           # C5: CH–OH
                   "[CX4H2][OX2H]")           # C6: CH2–OH
    keto_query = Chem.MolFromSmarts(keto_smarts)
    if keto_query is None:
        return None, None

    # Check open-chain aldo–hexose match.
    matches = mol.GetSubstructMatches(aldo_query)
    # We expect the SMARTS to have a number of atoms equal to the query’s atom count.
    if matches:
        # In our query, the total number of atoms is:
        n_query_atoms = aldo_query.GetNumAtoms()  # should be 12 (6 carbons + 6 oxygens)
        for match in matches:
            if len(match) == n_query_atoms:
                # Additionally, a rough check: count how many carbon atoms in the match.
                match_carbons = sum(1 for idx in match if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
                if match_carbons == 6:
                    return True, "Matches open‐chain aldo‐hexose pattern (aldehyde at C1)"
    
    # Check open–chain keto–hexose match.
    matches = mol.GetSubstructMatches(keto_query)
    if matches:
        n_query_atoms = keto_query.GetNumAtoms()  # should be 12 as well.
        for match in matches:
            if len(match) == n_query_atoms:
                match_carbons = sum(1 for idx in match if mol.GetAtomWithIdx(idx).GetSymbol() == "C")
                if match_carbons == 6:
                    return True, "Matches open‐chain keto‐hexose pattern (ketone at C2)"
    
    # -------------------------
    # 3. Cyclic (ring) detection.
    #
    # We look for either a pyranose–like ring (6–membered ring with exactly one oxygen, i.e. 5 C)
    # or a furanose–like ring (5–membered ring with exactly one oxygen, i.e. 4 C) that further bears an
    # exocyclic CH2OH branch so that the entire sugar core comprises 6 carbons.
    ring_info = mol.GetRingInfo()
    total_heavy = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    for ring in ring_info.AtomRings():
        # Look for pyranose candidate: ring of 6 atoms.
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            n_ox_ring = sum(1 for atom in ring_atoms if atom.GetSymbol() == "O")
            if n_ox_ring != 1:
                continue
            n_C_in_ring = sum(1 for atom in ring_atoms if atom.GetSymbol() == "C")
            # For a hexose in pyranose form, we expect 5 C atoms in the ring.
            if n_C_in_ring != 5:
                continue
            # Now look for an exocyclic CH2OH substituent attached to one of the ring carbons.
            exocyclic_found = False
            exo_idx = None
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() != "C":
                    continue
                # For each neighbor not in the ring…
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    # We want a CH2 group (approximately, at least 2 attached Hs) that is connected to an –OH.
                    if nbr.GetSymbol() == "C" and nbr.GetTotalNumHs() >= 2:
                        # Check if this carbon has an oxygen neighbor (at least one H on it)
                        for subnbr in nbr.GetNeighbors():
                            if subnbr.GetIdx() == atom.GetIdx():
                                continue
                            if subnbr.GetSymbol() == "O" and subnbr.GetTotalNumHs() >= 1:
                                exocyclic_found = True
                                exo_idx = nbr.GetIdx()
                                break
                    if exocyclic_found:
                        break
                if exocyclic_found:
                    break
            if exocyclic_found:
                # Now the candidate sugar core comprises the ring (5 C) plus the exocyclic CH2OH (1 C) = 6 C total.
                # As a safeguard, we require that the candidate substructure is not an insignificant decoration:
                # we require that the sum of atoms in the ring plus branch is at least 30% of the heavy atoms.
                candidate_heavy = len(ring) + 2  # roughly (1 extra C and its attached O)
                if candidate_heavy >= 0.3 * total_heavy:
                    return True, "Contains pyranose ring pattern (6–membered ring with one oxygen and an exocyclic CH2OH branch)"
        
        # Furanose candidate: ring of 5 atoms.
        if len(ring) == 5:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            n_ox_ring = sum(1 for atom in ring_atoms if atom.GetSymbol() == "O")
            if n_ox_ring != 1:
                continue
            n_C_in_ring = sum(1 for atom in ring_atoms if atom.GetSymbol() == "C")
            if n_C_in_ring != 4:
                continue
            # Look for an exocyclic CH2OH branch on one of the ring carbons.
            exocyclic_found = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() != "C":
                    continue
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    if nbr.GetSymbol() == "C" and nbr.GetTotalNumHs() >= 2:
                        for subnbr in nbr.GetNeighbors():
                            if subnbr.GetIdx() == atom.GetIdx():
                                continue
                            if subnbr.GetSymbol() == "O" and subnbr.GetTotalNumHs() >= 1:
                                exocyclic_found = True
                                break
                    if exocyclic_found:
                        break
                if exocyclic_found:
                    break
            if exocyclic_found:
                candidate_heavy = len(ring) + 2  # roughly candidate has 4+1 = 5 C? Actually for a genuine hexose furanose, it should total 5 C,
                                                # so we check that the branch brings the total to 5; but then a hexose must have 6 C.
                                                # Many hexoses occur in furanose form (e.g. tagatofuranose) where the CH2OH makes up the 6th C.
                # Here we assume that the furanose ring plus its exocyclic branch gives the hexose core.
                if candidate_heavy >= 0.3 * total_heavy:
                    return True, "Contains furanose ring pattern (5–membered ring with one oxygen and an exocyclic CH2OH branch)"
    
    # If no strategy succeeded:
    return False, "Does not match recognized hexose patterns"

# ----- Example usage -----
if __name__ == '__main__':
    # These test cases include some of the provided positive and problematic cases.
    test_smiles = [
        "OC(C(O)CNCCCCCCC)C(O)C(O)CO",                         # 1-Deoxy-1-(heptylamino)hexitol (modified hexose; open-chain variant)
        "O=C(OC1OC(C(O)C(C1O)O)C)C2=CC=CC=C2",                 # 1-O-Benzoyl-alpha-L-rhamnopyranoside (cyclic with acyl substitution)
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",               # D-allopyranose (cyclic)
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",         # aldehydo-D-galactose (open-chain)
        "C[C@H](O)[C@H](O)[C@@H](O)C(=O)CO",                   # L-rhamnulose (open-chain keto variant)
        "[H]C([H])([C@]([H])(O)C=O)[C@]([H])(O)[C@@]([H])(C)O",  # tyvelose (open-chain with deoxy substitution)
        "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H]1O",               # alpha-D-tagatofuranose (cyclic furanose)
    ]
    for s in test_smiles:
        result, reason = is_hexose(s)
        print(f"SMILES: {s}\n  Classified as hexose? {result}. Reason: {reason}\n")