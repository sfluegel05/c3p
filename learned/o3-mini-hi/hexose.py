"""
Classifies: CHEBI:18133 hexose
"""
#!/usr/bin/env python
"""
Classifies: hexose
Definition: Any six‐carbon monosaccharide which, in its linear form, contains either an aldehyde group
at position 1 (aldohexose) or a ketone group at position 2 (ketohexose). In its cyclic form the sugar
often forms a pyranose (6–membered ring: 5 C + 1 O plus an exocyclic CH2OH) or a furanose (5–membered ring
plus exocyclic CH2OH group(s)) ring. (All criteria are heuristic.)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule qualifies as a hexose.
    
    Strategy (heuristic):
      1. Quick rejection: if the molecule contains phosphorus.
      2. Open–chain detection:
         Look for an aldo–hexose pattern (roughly: O=CH-(CHOH)4-CH2OH) or a keto–hexose 
         (roughly: CH2OH-C(=O)-(CHOH)3-CH2OH). Then verify that exactly 6 C atoms are involved 
         and that extra substituents (attached as C–bonds) are minimal.
      3. Cyclic detection:
         For a pyranose candidate, look for a 6–membered ring that has 5 C and 1 O plus exactly one 
         exocyclic carbon (which should be –CH2OH). For a furanose candidate, require a 5–membered ring 
         (4 C and 1 O) plus two exocyclic carbons. In both cases we then check that apart from the expected 
         exocyclic carbon(s), the ring atoms have no extra carbon attachments.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as hexose, else False.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # add explicit hydrogens so we can count H atoms reliably    
    mol = Chem.AddHs(mol)
    
    # 1. Quick rejection: if molecule contains phosphorus.
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Contains phosphorus, likely not a hexose"

    # --- Helper functions ---
    def count_atoms(atom_indices, symbol):
        return sum(1 for idx in atom_indices if mol.GetAtomWithIdx(idx).GetSymbol() == symbol)
    
    def extra_carbon_penalty(atom_idx, core_set):
        """
        Count how many neighboring carbons (that are not in the core set) are attached.
        """
        count = 0
        for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            if nbr.GetSymbol() == "C" and nbr.GetIdx() not in core_set:
                count += 1
        return count
    
    def check_open_chain_match(match):
        """
        Given a substructure match (tuple of atom indices) for an open-chain pattern,
        verify that exactly 6 carbons are in the match and that no carbon in the match carries
        an extra carbon substituent.
        """
        core = set(match)
        carbons = [idx for idx in core if mol.GetAtomWithIdx(idx).GetSymbol() == "C"]
        if len(carbons) != 6:
            return False
        # For an ideal open-chain hexose, each carbon should only be attached (within this chain) to two others,
        # except for the ends (attached to one) and, possibly, the anomeric position. Extra carbon substituents (coming
        # from, e.g., long chains) are a red flag.
        penalty = 0
        for idx in carbons:
            penalty += extra_carbon_penalty(idx, core)
        # We allow a penalty of 0 or 1 (for common modifications, e.g. deoxy or O-glycosidic substitution at C1)
        return penalty <= 1

    def check_pyranose_match(ring):
        """
        Given a candidate ring (list of atom indices), check for a pyranose sugar core.
        Expect 6 member atoms: 5 carbons and 1 oxygen.
        Then look for exactly one exocyclic carbon attached to a ring carbon.
        Also require that no other ring carbon has an additional carbon neighbor.
        """
        ring_set = set(ring)
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        if symbols.count("C") != 5 or symbols.count("O") != 1:
            return False
        # Now collect exocyclic carbon neighbors attached to ring carbons.
        exo = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C": 
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetSymbol() == "C":
                    exo.add(nbr.GetIdx())
        if len(exo) != 1:
            return False
        # The total number of carbons in the putative sugar is 5 (in ring) + 1 (exo) = 6.
        # Check that the exocyclic carbon appears to be a CH2OH group:
        exo_atom = mol.GetAtomWithIdx(list(exo)[0])
        if exo_atom.GetTotalNumHs() < 2:
            return False
        if not any(nbr.GetSymbol() == "O" for nbr in exo_atom.GetNeighbors() if nbr.GetIdx() not in ring_set):
            return False
        # Finally, check that none of the ring carbons (apart from the one bearing the exocyclic CH2OH)
        # has any extra carbon attachments.
        exo_attached = list(exo)[0]
        allowed = {exo_attached}
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C": 
                continue
            extras = [nbr.GetIdx() for nbr in atom.GetNeighbors() 
                      if nbr.GetSymbol() == "C" and nbr.GetIdx() not in ring_set]
            # Allow the one expected exocyclic neighbor only once.
            if len(extras) > (1 if allowed.intersection(extras) else 0):
                return False
        return True

    def check_furanose_match(ring):
        """
        For a furanose candidate we expect a 5–member ring with 4 carbons and 1 oxygen,
        plus two exocyclic carbons (each expected to be CH2OH) so that total C count is 6.
        Also allow no extra C attachments.
        """
        ring_set = set(ring)
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        if symbols.count("C") != 4 or symbols.count("O") != 1:
            return False
        exo = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_set:
                    continue
                if nbr.GetSymbol() == "C":
                    exo.add(nbr.GetIdx())
        if len(exo) != 2:
            return False
        # Check each exocyclic carbon looks like CH2OH.
        for exo_idx in exo:
            ea = mol.GetAtomWithIdx(exo_idx)
            if ea.GetTotalNumHs() < 2:
                return False
            if not any(nbr.GetSymbol() == "O" for nbr in ea.GetNeighbors() if nbr.GetIdx() not in ring_set):
                return False
        # Check that none of the ring carbons (other than those bonded to the exo groups) has extra carbon attachments.
        allowed = exo
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            extras = [nbr.GetIdx() for nbr in atom.GetNeighbors() 
                      if nbr.GetSymbol() == "C" and nbr.GetIdx() not in ring_set]
            # each ring carbon should only have the exo substituents that we expect (and they may be different on different atoms)
            if len(extras) > (1 if set(extras).intersection(allowed) else 0):
                return False
        return True

    # --- 2. Open–chain detection using relaxed SMARTS.
    # Aldo–hexose: roughly O=CH-(CHOH)4-CH2OH (relaxed to allow substitutions)
    aldo_smarts = ("[CX3H1](=O)"       # aldehyde carbon (it should have one hydrogen)
                   "[CX4]"             # then a CH (bearing an -OH) 
                   "[CX4]"             # repeated three times
                   "[CX4]"
                   "[CX4]"             # then terminal CH2OH (hidden: we allow it to be C without explicit H count)
                   )
    aldo_query = Chem.MolFromSmarts(aldo_smarts)
    if aldo_query is None:
        return None, None  # Should not happen
    for match in mol.GetSubstructMatches(aldo_query):
        if check_open_chain_match(match):
            return True, "Matches open–chain aldo–hexose pattern (aldehyde at C1)"
    
    # Keto–hexose: roughly CH2OH-C(=O)-(CHOH)3-CH2OH
    keto_smarts = ("[CX4]"             # terminal CH2OH at one end
                   "[CX3](=O)"         # ketone carbon
                   "[CX4]"             # then CHOH's
                   "[CX4]"
                   "[CX4]"
                   "[CX4]"             # terminal CH2OH at the other end
                   )
    keto_query = Chem.MolFromSmarts(keto_smarts)
    if keto_query is None:
        return None, None
    for match in mol.GetSubstructMatches(keto_query):
        if check_open_chain_match(match):
            return True, "Matches open–chain keto–hexose pattern (ketone at C2)"
    
    # --- 3. Cyclic (ring) detection.
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    for ring in rings:
        if len(ring) == 6:  # pyranose candidate: expect 5 C + 1 O in ring.
            if check_pyranose_match(ring):
                return True, "Contains pyranose ring pattern (6-membered ring with one oxygen and one exocyclic CH2OH branch)"
        if len(ring) == 5:  # furanose candidate: expect 4 C + 1 O in ring plus 2 exocyclic C's.
            if check_furanose_match(ring):
                return True, "Contains furanose ring pattern (5-membered ring with one oxygen and sufficient exocyclic carbons)"
    
    return False, "Does not match recognized hexose patterns"

# ----- Example usage -----
if __name__ == '__main__':
    test_smiles = [
        "OC(C(O)CNCCCCCCC)C(O)C(O)CO",                         # 1-Deoxy-1-(heptylamino)hexitol
        "O=C(OC1OC(C(O)C(C1O)O)C)C2=CC=CC=C2",                 # 1-O-Benzoyl-alpha-L-rhamnopyranoside
        "OC[C@H]1OC(O)[C@H](O)[C@H](O)[C@@H]1O",               # D-allopyranose
        "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)CO",         # aldehydo-D-galactose
        "C[C@H](O)[C@H](O)[C@@H](O)C(=O)CO",                   # L-rhamnulose
        "[H][C@]([H])(O)C=O[C@]([H])(O)[C@@]([H])(C)O",         # tyvelose (open-chain with deoxy substitution)
        "OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H]1O",               # alpha-D-tagatofuranose
    ]
    for s in test_smiles:
        res, reason = is_hexose(s)
        print(f"SMILES: {s}\n  Classified as hexose? {res}. Reason: {reason}\n")