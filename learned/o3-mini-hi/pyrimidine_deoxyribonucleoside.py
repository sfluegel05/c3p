"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
#!/usr/bin/env python
"""
Classifies: pyrimidine deoxyribonucleoside

A pyrimidine deoxyribonucleoside is defined as a deoxyribonucleoside containing a pyrimidine base.
That is, it comprises a deoxyribose sugar (typically a five‐membered (furanose) ring with four carbons and one oxygen,
with one fewer exocyclic –OH than ribose) connected via a glycosidic C–N bond to a pyrimidine base.
The pyrimidine base (a six‐membered ring) should have at least two nitrogen atoms and exhibit a carbonyl (C=O).
Molecules containing phosphorus are rejected.
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines whether a molecule is a pyrimidine deoxyribonucleoside 
    based on its SMILES string.
    
    The method attempts to identify a candidate deoxyribose sugar – a 5‐membered (furanose) ring 
    made of exactly 4 carbons and 1 oxygen that has only 1–2 exocyclic oxygen substituents 
    (the lower count compared with ribose, which normally carries 3 such substituents) – 
    and a candidate pyrimidine base – a 6‐membered ring with at least 2 nitrogens and a carbonyl group.
    It then verifies that these two substructures are connected via a single glycosidic (C–N) bond.
    
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
    
    # Reject any molecules with phosphorus (e.g. nucleotides).
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus and is likely not a deoxyribonucleoside."
    
    # Add explicit hydrogens (needed for hydroxyl group detection).
    mol = Chem.AddHs(mol)
    
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Each ring is a tuple of atom indices.
    
    # STEP 1: Detect candidate deoxyribose sugar ring.
    # We require a 5-membered ring (furanose) made up of exactly 4 carbons and 1 oxygen.
    # In 2'-deoxyribose (or in some dideoxy sugars) the count of exocyclic oxygen substituents (attached to ring carbons)
    # is lower (typically 2 or even 1) compared with ribose (3 such oxygens).
    sugar_candidate = None
    sugar_reason = "No candidate deoxyribose sugar ring found."
    for ring in atom_rings:
        if len(ring) != 5:
            continue
        # Count number of carbon and oxygen atoms in the ring.
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
        if not valid_ring or c_count != 4 or o_count != 1:
            continue
        
        # Count the exocyclic oxygen substituents (as candidates for –OH or –CH2OH).
        # (Only consider substituents attached to ring carbons.)
        exo_o_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                # Only consider neighbors outside the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Look for oxygen atoms that are likely to be part of an -OH or -CH2OH.
                if nbr.GetSymbol() == "O":
                    # Check that the oxygen carries at least one hydrogen.
                    if any(neigh.GetSymbol() == "H" for neigh in nbr.GetNeighbors()):
                        exo_o_count += 1
        # For a deoxyribose (or dideoxy sugar) we expect an exo_oxygen count of 1 or 2,
        # whereas ribose should have 3.
        if exo_o_count in (1, 2):
            sugar_candidate = set(ring)
            sugar_reason = (f"Found candidate deoxyribose sugar ring: 5-membered ring with 4 C + 1 O "
                            f"and {exo_o_count} exocyclic oxygen substituent(s) (lower than ribose).")
            break
    if sugar_candidate is None:
        return False, sugar_reason
    
    # STEP 2: Detect candidate pyrimidine base.
    # We require a 6-membered ring that contains at least 2 nitrogen atoms and evidence of a carbonyl group.
    base_candidate = None
    base_reason = "No candidate pyrimidine base (6-membered ring with ≥2 N and a carbonyl group) found."
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        n_count = sum(1 for atom in atoms_in_ring if atom.GetSymbol() == "N")
        if n_count < 2:
            continue
        
        # Look for a carbonyl group.
        carbonyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            # Check bonds from this carbon. A double bond to oxygen is taken as a carbonyl.
            for bond in atom.GetBonds():
                if bond.GetBondTypeAsDouble() == 2.0:
                    # Identify the neighboring atom.
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

    # STEP 3: Ensure the sugar and base are connected by a single glycosidic bond.
    # In a (deoxy)ribonucleoside the connection is typically a C–N bond between a sugar carbon and a base nitrogen.
    glyco_connections = []
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        if ((idx1 in sugar_candidate and idx2 in base_candidate) or 
            (idx2 in sugar_candidate and idx1 in base_candidate)):
            atom1 = mol.GetAtomWithIdx(idx1)
            atom2 = mol.GetAtomWithIdx(idx2)
            if ((atom1.GetSymbol() == "C" and atom2.GetSymbol() == "N") or 
                (atom1.GetSymbol() == "N" and atom2.GetSymbol() == "C")):
                glyco_connections.append((idx1, idx2))
    if len(glyco_connections) != 1:
        return False, ("Sugar and base are not connected by a single glycosidic C–N bond "
                       f"(found {len(glyco_connections)} such connection(s)).")
    
    # If all criteria are met, return True with an explanation.
    return True, ("Molecule contains a deoxyribose sugar (5-membered ring with 4 C + 1 O and 1–2 exocyclic –OH/–CH2OH substituents) "
                  "connected via a glycosidic C–N bond to a pyrimidine base (6-membered ring with ≥2 N and a carbonyl group).")

# (Optional) Quick testing if the script is run directly.
if __name__ == '__main__':
    # A small selection of examples: true positives and an example that was previously missed.
    test_examples = [
        # True positives:
        ("S(N1C(=O)N(C=C(C1=O)C)[C@@H]2O[C@H](CO)[C@H](C2)O)N", "Thymidine-3-thioamine"),
        ("O1[C@@H](N2C(O)C(O)(C(=O)NC2=O)C)C[C@H](O)[C@H]1CO", "Thymidine glycol"),
        ("[C@@H]1(N2C(NC(=O)C(=C2)C(=O)O)=O)O[C@H](CO)[C@H](C1)O", "5-carboxy-2'-deoxyuridine"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)nc1N", "5-methyl-2'-deoxycytidine"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O", "thymidine"),
        # Previously missed example (false negative):
        ("FC1=CN(C2OC(CC2)CO)C(=O)N=C1N", "Endalin"),
        # Example that in previous attempt was missed:
        ("Nc1ccn([C@H]2CC[C@@H](CO)O2)c(=O)n1", "zalcitabine"),
    ]
    
    for smi, name in test_examples:
        res, reason = is_pyrimidine_deoxyribonucleoside(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*60}")