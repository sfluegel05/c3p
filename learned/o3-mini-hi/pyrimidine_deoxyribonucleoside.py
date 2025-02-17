"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
#!/usr/bin/env python
"""
Classifies: pyrimidine deoxyribonucleoside

A pyrimidine deoxyribonucleoside is defined as a deoxyribonucleoside containing a pyrimidine base.
In our approach the molecule must contain a candidate deoxyribose sugar – typically a 5‐membered (or sometimes 4‐membered)
ring containing one oxygen and four (or three) carbons with a low count of exocyclic –OH/–CH2OH substituents –
connected by a glycosidic C–N bond to a pyrimidine base (a 6‐membered ring with at least two nitrogens and a carbonyl group).
Molecules containing phosphorus are rejected.
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines whether a molecule is a pyrimidine deoxyribonucleoside 
    based on its SMILES string.

    The method searches for a candidate deoxyribose sugar – meaning a 5‐membered ring (or in some cases 4‐membered)
    containing exactly one oxygen (and respectively 4 or 3 carbons) that has a low number of exocyclic oxygen substituents,
    and a candidate pyrimidine base – a 6‐membered ring with at least 2 nitrogen atoms and a carbonyl (C=O)
    group. It finally requires that the two substructures be connected by a single C–N (glycosidic) bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a pyrimidine deoxyribonucleoside, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES and check for validity.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Reject molecules with phosphorus (e.g. nucleotides/phosphates)
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus and is likely not a deoxyribonucleoside."
    
    # Add explicit hydrogens (helpful for identifying hydroxyl groups).
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # STEP 1: Find candidate deoxyribose sugar ring.
    # We allow either a 5-membered ring (expect 4 carbons + 1 oxygen) or a 4-membered ring (3 carbons + 1 oxygen)
    sugar_candidate = None
    sugar_reason = "No candidate deoxyribose (or deoxy‐like) sugar ring found."
    for ring in atom_rings:
        n = len(ring)
        if n not in (5, 4):
            continue
        # Count atoms in the ring.
        c_count = 0
        o_count = 0
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            sym = atom.GetSymbol()
            if sym == "C":
                c_count += 1
            elif sym == "O":
                o_count += 1
            else:
                valid_ring = False
                break
        if not valid_ring:
            continue
        # For a 5-membered ring we expect 4 Cs, 1 O; for a 4-membered ring we expect 3 Cs, 1 O.
        if (n == 5 and (c_count, o_count) != (4, 1)) or (n == 4 and (c_count, o_count) != (3, 1)):
            continue

        # Count exocyclic oxygen substituents on ring carbons.
        exo_o_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms in the ring.
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "O":
                    # Only count if oxygen is not part of a carbonyl (i.e. it has a hydrogen attached)
                    if any(nh.GetSymbol() == "H" for nh in nbr.GetNeighbors()):
                        exo_o_count += 1
        # For a deoxyribose sugar (5-membered) we expect typically 1 or 2 exocyclic O; 
        # for a four-membered sugar analog, we allow a slightly lower count (0 to 2).
        if (n == 5 and exo_o_count in (1, 2)) or (n == 4 and exo_o_count in (0, 1, 2)):
            sugar_candidate = set(ring)
            sugar_reason = (f"Found candidate sugar ring of size {n} with {c_count} C and {o_count} O, "
                            f"and {exo_o_count} exocyclic oxygen substituent(s).")
            break
    if sugar_candidate is None:
        return False, sugar_reason

    # STEP 2: Find candidate pyrimidine base.
    # Look for a 6-membered ring with at least 2 nitrogen atoms
    # and with evidence of a carbonyl (C=O) group.
    base_candidate = None
    base_reason = "No candidate pyrimidine base (6-membered ring with ≥2 N and at least one carbonyl) found."
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        n_nitrogen = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == "N")
        if n_nitrogen < 2:
            continue
        carbonyl_found = False
        # For each carbon in the ring see if it is double-bonded to an oxygen.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for bond in atom.GetBonds():
                # Check if the bond is a double bond.
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetSymbol() == "O":
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        if not carbonyl_found:
            continue
        
        base_candidate = set(ring)
        base_reason = "Found candidate pyrimidine base: 6-membered ring with ≥2 N and a carbonyl group."
        break
    if base_candidate is None:
        return False, base_reason

    # STEP 3: Check that the sugar and base are connected by a single glycosidic C–N bond.
    glyco_connections = []
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        # Connection must be between an atom from sugar and an atom from base.
        if ((idx1 in sugar_candidate and idx2 in base_candidate) or 
            (idx2 in sugar_candidate and idx1 in base_candidate)):
            # And the bond should be between C and N.
            a1 = mol.GetAtomWithIdx(idx1)
            a2 = mol.GetAtomWithIdx(idx2)
            if ((a1.GetSymbol() == "C" and a2.GetSymbol() == "N") or 
                (a1.GetSymbol() == "N" and a2.GetSymbol() == "C")):
                glyco_connections.append((idx1, idx2))
    if len(glyco_connections) != 1:
        return False, ("Sugar and base are not connected by a single glycosidic C–N bond "
                       f"(found {len(glyco_connections)} connection(s)).")
    
    return True, ("Molecule contains a deoxyribose (or deoxy‐like) sugar ring "
                  f"({sugar_reason}) connected via a glycosidic C–N bond to a pyrimidine base "
                  f"({base_reason}).")

# (Optional) Quick testing if run directly.
if __name__ == '__main__':
    test_examples = [
        # True positives:
        ("S(N1C(=O)N(C=C(C1=O)C)[C@@H]2O[C@H](CO)[C@H](C2)O)N", "Thymidine-3-thioamine"),
        ("O1[C@@H](N2C(O)C(O)(C(=O)NC2=O)C)C[C@H](O)[C@H]1CO", "Thymidine glycol"),
        ("[C@@H]1(N2C(NC(=O)C(=C2)C(=O)O)=O)O[C@H](CO)[C@H](C1)O", "5-carboxy-2'-deoxyuridine"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)nc1N", "5-methyl-2'-deoxycytidine"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O", "thymidine"),
        ("N1(C(NC(C(=C1)CCN)=O)=O)[C@H]2C[C@@H]([C@H](O2)CO)O", "5-(2-aminoethyl)-2'-deoxyuridine"),
        ("OC[C@H]1O[C@H](C[C@@H]1O)n1cc(Br)c(=O)[nH]c1=O", "5-bromo-2'-deoxyuridine"),
        ("[3H]Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O", "alpha-tritiated thymidine"),
        # Previously missed examples:
        ("FC1=CN(C2OC(CC2)CO)C(=O)N=C1N", "Endalin"),
        ("Nc1ccn([C@H]2CC[C@@H](CO)O2)c(=O)n1", "zalcitabine"),
    ]
    for smi, name in test_examples:
        res, reason = is_pyrimidine_deoxyribonucleoside(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")