"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
#!/usr/bin/env python
"""
Classifies: pyrimidine deoxyribonucleoside

A pyrimidine deoxyribonucleoside is defined as a deoxyribonucleoside containing a pyrimidine base.
That is, it comprises a deoxyribose sugar (a five‐membered ring made of four carbons and one oxygen, 
with the characteristic absence of a 2′-OH – meaning that only one hydroxyl group is directly present on the ring)
connected via a glycosidic C–N bond to a pyrimidine base. The base candidate is a six–membered ring 
which should contain at least two nitrogen atoms and evidence of a carbonyl (C=O) group. Molecules 
containing phosphorus are rejected.
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
    
    # Reject molecules containing phosphorus (e.g. nucleotides).
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus atoms and is likely not a deoxyribonucleoside."
    
    # Add explicit hydrogens to assist in detecting hydroxyl groups.
    mol = Chem.AddHs(mol)
    
    # Obtain ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    sugar_candidate = None
    sugar_reason = "No candidate deoxyribose sugar ring found."
    
    # STEP 1: Identify candidate deoxyribose sugar rings.
    # We look for 5-membered rings (furanose) that contain exactly 4 carbon atoms and 1 oxygen.
    # Furthermore, among the ring carbons, we require exactly one exocyclic hydroxyl (an O not in the ring
    # that is attached to a ring carbon and has at least one H neighbor). This reflects the absence of a 2'-OH.
    for ring in atom_rings:
        if len(ring) != 5:
            continue
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
        
        # Count exocyclic hydroxyl groups on ring carbon atoms.
        hydroxy_exocyclic = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            # Look at neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # If the neighbor is oxygen and is attached to at least one hydrogen then count it.
                if nbr.GetSymbol() == "O":
                    if any(n.GetSymbol() == "H" for n in nbr.GetNeighbors()):
                        hydroxy_exocyclic += 1
        # For a deoxyribose, we expect exactly one oxygen directly attached as an -OH on a ring carbon.
        if hydroxy_exocyclic == 1:
            sugar_candidate = set(ring)
            sugar_reason = ("Found candidate deoxyribose sugar ring: 5-membered ring with 4 C + 1 O "
                            "and exactly 1 exocyclic hydroxyl (reflecting absence of 2'-OH).")
            break
    if sugar_candidate is None:
        return False, sugar_reason
    
    # STEP 2: Identify candidate pyrimidine base rings.
    # Look for 6-membered rings with at least 2 nitrogen atoms and evidence of a carbonyl group.
    base_candidate = None
    base_reason = "No candidate pyrimidine base (6-membered ring with ≥2 N and a carbonyl group) found."
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "N")
        if n_count < 2:
            continue
        
        # Look for a carbonyl group either as a double bond within the ring or as an exocyclic substituent.
        carbonyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for bond in atom.GetBonds():
                # If bond order is double (i.e. a double bond) and the neighboring atom is oxygen.
                if bond.GetBondTypeAsDouble() == 2.0:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetSymbol() == "O":
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        if not carbonyl_found:
            continue
        
        base_candidate = set(ring)
        base_reason = ("Found candidate pyrimidine base: 6-membered ring with at least 2 N and a carbonyl group.")
        break
    if base_candidate is None:
        return False, base_reason
    
    # STEP 3: Ensure the sugar and base are connected by a glycosidic bond.
    # In a deoxyribonucleoside, the connection is typically a C–N bond between a sugar carbon and a base nitrogen.
    glyco_bond_found = False
    for bond in mol.GetBonds():
        a1_idx = bond.GetBeginAtomIdx()
        a2_idx = bond.GetEndAtomIdx()
        if ((a1_idx in sugar_candidate and a2_idx in base_candidate) or
            (a2_idx in sugar_candidate and a1_idx in base_candidate)):
            a1 = mol.GetAtomWithIdx(a1_idx)
            a2 = mol.GetAtomWithIdx(a2_idx)
            if (a1.GetSymbol() == "C" and a2.GetSymbol() == "N") or (a1.GetSymbol() == "N" and a2.GetSymbol() == "C"):
                glyco_bond_found = True
                break
    if not glyco_bond_found:
        return False, "Sugar and base are not connected via a glycosidic (C–N) bond."
    
    return True, ("Molecule contains a deoxyribose sugar (5-membered ring of 4 C + 1 O with a single exocyclic –OH) "
                  "joined via a glycosidic C–N bond to a pyrimidine base (6-membered ring with ≥2 N and a carbonyl group).")

# (Optional) Quick testing when run as a script.
if __name__ == '__main__':
    test_examples = [
        # True positives
        ("S(N1C(=O)N(C=C(C1=O)C)[C@@H]2O[C@H](CO)[C@H](C2)O)N", "Thymidine-3-thioamine"),
        ("O1[C@@H](N2C(O)C(O)(C(=O)NC2=O)C)C[C@H](O)[C@H]1CO", "Thymidine glycol"),
        ("[C@@H]1(N2C(NC(=O)C(=C2)C(=O)O)=O)O[C@H](CO)[C@H](C1)O", "5-carboxy-2'-deoxyuridine"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)nc1N", "5-methyl-2'-deoxycytidine"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O", "thymidine"),
        # Examples that were false negatives before:
        ("FC1=CN(C2OC(CC2)CO)C(=O)N=C1N", "Endalin"),
    ]
    for smi, name in test_examples:
        res, reason = is_pyrimidine_deoxyribonucleoside(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")