"""
Classifies: CHEBI:17334 penicillin
"""
#!/usr/bin/env python
"""
Classifies: Penicillins
Definition:
  “Any member of the group of substituted penams containing two methyl substituents at position 2,
   a carboxylate substituent at position 3 and a carboxamido group at position 6.”
   
This implementation first removes potential counterions by selecting the largest fragment,
detects the fused bicyclic penicillin core, and then restricts subgroup searches (for the dimethyl,
carboxylate, and carboxamido groups) to atoms that are part of (or directly attached to) the detected core.
"""

from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    
    A molecule is considered a penicillin if:
      (a) Its largest organic fragment (after salt removal) is neutral.
      (b) It contains a fused bicyclic penicillin core – namely, a 4-membered β‐lactam (with a nitrogen
          and at least one carbonyl) fused to a 5-membered thiazolidine (containing at least one sulfur).
      (c) In that core the carbon adjancent to the sulfur (position 2) bears two methyl substituents.
      (d) A carboxylate substituent (C(=O)O or C(=O)[O-]) is attached to the core.
      (e) A carboxamido fragment (N–C(=O)) is attached to the core.
      
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule appears to be a penicillin; False otherwise.
        str: Explanation of the decision.
    """
    # Parse molecule from SMILES and isolate the largest fragment (to remove counterions/salts)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if not frags:
        return False, "No fragments found in molecule"
    frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    if Chem.GetFormalCharge(frag) != 0:
        return False, "Major fragment is charged (likely a salt) and does not count as a penicillin core"

    # --- Step 1: Identify the fused penicillin core ---
    # We require:
    #   (i) A 4-membered ring (β-lactam) containing at least one N and one C that is double bonded to O.
    #   (ii) A 5-membered ring (thiazolidine) that shares exactly two atoms with the 4-membered ring and has at least one S.
    def get_penicillin_core(mol):
        ring_info = mol.GetRingInfo()
        rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
        for ring4 in rings:
            if len(ring4) != 4:
                continue
            ring4_set = set(ring4)
            atoms4 = [mol.GetAtomWithIdx(i) for i in ring4_set]
            if not any(atom.GetSymbol() == "N" for atom in atoms4):
                continue
            has_carbonyl = False
            for atom in atoms4:
                if atom.GetSymbol() == "C":
                    for bond in atom.GetBonds():
                        if bond.GetBondType().name == "DOUBLE":
                            neigh = bond.GetOtherAtom(atom)
                            if neigh.GetSymbol() == "O":
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        break
            if not has_carbonyl:
                continue
            # Look for a 5-membered ring sharing exactly 2 atoms with the 4-membered ring and containing S
            for ring5 in rings:
                if len(ring5) != 5:
                    continue
                common = ring4_set.intersection(ring5)
                if len(common) == 2:
                    atoms5 = [mol.GetAtomWithIdx(i) for i in ring5]
                    if any(atom.GetSymbol() == "S" for atom in atoms5):
                        return ring4_set.union(ring5)
        return None

    core_atoms = get_penicillin_core(frag)
    if core_atoms is None:
        return False, "Molecule does not contain the required fused penicillin core (β-lactam fused with thiazolidine containing sulfur)"

    # --- Step 2: Check for two methyl substituents on the carbon adjacent to the core sulfur (position 2) ---
    # For any sulfur within the core, search its neighboring carbon in the core and count methyl groups (CH3) attached,
    # ensuring they lie outside of the core.
    found_dimethyl = False
    for atom in frag.GetAtoms():
        if atom.GetIdx() not in core_atoms:
            continue
        if atom.GetSymbol() != "S":
            continue
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in core_atoms or nbr.GetSymbol() != "C":
                continue
            methyl_count = 0
            for sub in nbr.GetNeighbors():
                if sub.GetIdx() in core_atoms:
                    continue
                if sub.GetAtomicNum() == 6 and sub.GetDegree() == 1:
                    methyl_count += 1
            if methyl_count == 2:
                found_dimethyl = True
                break
        if found_dimethyl:
            break
    if not found_dimethyl:
        return False, "Missing two methyl substituents on the core carbon adjacent to the sulfur (position 2)"
    
    # --- Step 3: Check for a carboxylate substituent attached to the core ---
    # Instead of forcing the carboxylate carbon to be a ring atom, we define a more flexible SMARTS:
    # "[CX3](=O)[O,OX1-]" which can match both cyclic and exocyclic carboxylate groups.
    carboxylate_smarts = Chem.MolFromSmarts("[CX3](=O)[O,OX1-]")
    if carboxylate_smarts is None:
        return False, "Error in carboxylate SMARTS"
    carboxylate_found = False
    for match in frag.GetSubstructMatches(carboxylate_smarts):
        # match[0] is the carboxylate carbon; check if it is directly attached to a core atom.
        c_idx = match[0]
        carboxylate_atom = frag.GetAtomWithIdx(c_idx)
        for nbr in carboxylate_atom.GetNeighbors():
            if nbr.GetIdx() in core_atoms:
                carboxylate_found = True
                break
        if carboxylate_found:
            break
    if not carboxylate_found:
        return False, "Missing a carboxylate substituent (C(=O)O or C(=O)[O-]) attached to the penicillin core"

    # --- Step 4: Check for a carboxamido fragment (N–C(=O)) attached to the core ---
    # We use a SMARTS for an amide moiety. The standard fragment "[N;R]-C(=O)" is used, and we require that one atom in the fragment is part of the core.
    carboxamido_smarts = Chem.MolFromSmarts("[N]-C(=O)")
    if carboxamido_smarts is None:
        return False, "Error in carboxamido SMARTS"
    carboxamido_found = False
    for match in frag.GetSubstructMatches(carboxamido_smarts):
        if any(idx in core_atoms for idx in match):
            carboxamido_found = True
            break
    if not carboxamido_found:
        return False, "Missing a carboxamido group (N-C(=O)) attached to the penicillin core"
    
    return True, "Molecule has a fused penicillin core with dimethyl at position 2, a carboxylate at position 3, and a carboxamido at position 6"

# Example usage:
if __name__ == "__main__":
    # Test with Penicillin K SMILES (as one example)
    test_smiles = "CCCCCCCC(=O)N[C@H]1[C@H]2SC(C)(C)[C@@H](N2C1=O)C(O)=O"
    result, reason = is_penicillin(test_smiles)
    print("Is penicillin?", result)
    print("Reason:", reason)