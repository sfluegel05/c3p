"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
#!/usr/bin/env python
"""
Classifies: pyrimidine deoxyribonucleoside

A pyrimidine deoxyribonucleoside is defined as a deoxyribonucleoside containing a pyrimidine base.
That is, it comprises a deoxyribose sugar (a five‐membered ring made of four carbons and one oxygen,
with the characteristic absence of a 2′-OH, so that only one –OH appears on the ring)
connected via a glycosidic bond (i.e. C–N bond) to a pyrimidine base.
We require that the base candidate is a six–membered ring which contains at least two nitrogen atoms
and has evidence for at least one carbonyl group (C=O) either on a ring atom or as exocyclic substituent.
Molecules containing phosphorus are rejected (to avoid nucleotides).
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a pyrimidine deoxyribonucleoside, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Reject molecules containing phosphorus (likely nucleotides or larger P-containing compounds).
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus atoms and is likely not a deoxyribonucleoside."
    
    # Add explicit hydrogens to help detecting hydroxyl groups.
    mol = Chem.AddHs(mol)
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    sugar_candidate = None
    sugar_reason = "No candidate deoxyribose sugar ring found."
    # STEP 1: Identify candidate deoxyribose sugar rings:
    # Look for five-membered rings (furanose) made of exactly 4 carbons and 1 oxygen.
    # Additionally, we require that one of the ring carbons has an exocyclic –OH (as expected in deoxyribose).
    for ring in atom_rings:
        if len(ring) != 5:
            continue
        # Count atoms in the ring.
        c_count = 0
        o_count = 0
        valid_ring = True
        for idx in ring:
            sym = mol.GetAtomWithIdx(idx).GetSymbol()
            if sym == "C":
                c_count += 1
            elif sym == "O":
                o_count += 1
            else:
                valid_ring = False
                break
        if not valid_ring or c_count != 4 or o_count != 1:
            continue
        
        # Check for hydroxyl group(s) on ring carbons.
        # We only consider neighbors that are not in the ring and must be an oxygen bound to at least one hydrogen.
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "O":
                    # Check if the oxygen is part of a hydroxyl (has at least one H neighbor)
                    if any(n.GetSymbol() == "H" for n in nbr.GetNeighbors()):
                        hydroxyl_count += 1
        # For deoxyribose, one expects exactly one –OH to be directly bound on the ring.
        if hydroxyl_count == 1:
            sugar_candidate = set(ring)
            sugar_reason = "Found candidate deoxyribose sugar ring (5-membered with 4 C + 1 O and 1 exocyclic hydroxyl)."
            break
    if sugar_candidate is None:
        return False, sugar_reason
    
    # STEP 2: Identify candidate pyrimidine base rings.
    # We search for rings with 6 atoms that contain at least 2 nitrogen atoms.
    # Rather than requiring full aromaticity, we allow reduced rings.
    base_candidate = None
    base_reason = "No candidate pyrimidine base (6-membered ring with ≥2 N and a carbonyl group) found."
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        # Count nitrogen atoms in the ring.
        n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "N")
        if n_count < 2:
            continue
        
        # Look for evidence of a carbonyl group:
        # Check each atom in the ring (if it is carbon) and see if it is double-bonded to an oxygen atom
        # that lies outside the ring (exocyclic) or is part of the ring.
        carbonyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for bond in atom.GetBonds():
                # Check for double bond to oxygen.
                if bond.GetBondTypeAsDouble() == 2.0:
                    nbr = bond.GetOtherAtom(atom)
                    # Accept if oxygen is either in the ring or exocyclic.
                    if nbr.GetSymbol() == "O":
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        if not carbonyl_found:
            continue
        
        base_candidate = set(ring)
        base_reason = "Found candidate pyrimidine base (6-membered ring with ≥2 N and a carbonyl group)."
        break
    if base_candidate is None:
        return False, base_reason
    
    # STEP 3: Check for a glycosidic bond connecting the sugar and the base.
    # In a deoxyribonucleoside the glycosidic bond is usually between a carbon in the sugar candidate and a nitrogen in the base candidate.
    glyco_bond_found = False
    for bond in mol.GetBonds():
        a1_idx = bond.GetBeginAtomIdx()
        a2_idx = bond.GetEndAtomIdx()
        if ((a1_idx in sugar_candidate and a2_idx in base_candidate) or
            (a2_idx in sugar_candidate and a1_idx in base_candidate)):
            a1 = mol.GetAtomWithIdx(a1_idx)
            a2 = mol.GetAtomWithIdx(a2_idx)
            # Require that the bond is between a carbon (from the sugar) and a nitrogen (from the base)
            if (a1.GetSymbol() == "C" and a2.GetSymbol() == "N") or (a1.GetSymbol() == "N" and a2.GetSymbol() == "C"):
                glyco_bond_found = True
                break
    if not glyco_bond_found:
        return False, "Sugar and base are not connected by a glycosidic (C–N) bond."
    
    return True, "Molecule contains a deoxyribose sugar (5-membered ring with 4 C + 1 O and 1 hydroxyl) connected via a glycosidic C–N bond to a pyrimidine base (6-membered ring with ≥2 N and carbonyl)."

# (Optional) Quick testing when run as a script.
if __name__ == '__main__':
    test_examples = [
        # True positives (e.g. thymidine, 5-methyl-2'-deoxycytidine)
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O", "thymidine"),
        ("OC[C@H]1O[C@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O", "2'-deoxyuridine"),
        ("FC1=CN(C2OC(CC2)CO)C(=O)N=C1N", "Endalin"),
        ("O1[C@@H](N2C(O)C(O)(C(=O)NC2=O)C)C[C@H](O)[C@H]1CO", "Thymidine glycol")
    ]
    for smi, name in test_examples:
        res, reason = is_pyrimidine_deoxyribonucleoside(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")