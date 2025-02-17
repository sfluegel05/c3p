"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
#!/usr/bin/env python
"""
Classifies: pyrimidine deoxyribonucleoside
A deoxyribonucleoside is defined as a nucleoside containing a deoxyribose (i.e., lacking a 2'-OH)
sugar linked via a glycosidic bond to a pyrimidine base (aromatic six‐membered ring with 2 nitrogens).
Examples include thymidine, 2'-deoxyuridine, and 2'-deoxycytidine.
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    
    This function uses a set of heuristic rules:
      1. The molecule must be parsed successfully.
      2. The molecule must contain a candidate sugar ring that is five-membered (furanose),
         containing exactly one oxygen atom and four carbon atoms.
      3. The candidate sugar ring must display a pattern consistent with deoxyribose.
         In ribose nucleosides three hydroxyl substituents on the sugar are observed (attached to ring carbons),
         whereas in deoxyribose only two such hydroxyl groups are present.
      4. The molecule must contain a pyrimidine base—a six-membered aromatic ring containing exactly two nitrogen atoms.
      5. There must be at least one bond connecting the sugar ring and the pyrimidine base (i.e. the glycosidic bond).
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule passes all checks, False otherwise.
        str: Reason for the classification decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to allow proper detection of –OH groups.
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # STEP 1: Find candidate deoxyribose sugar rings.
    sugar_candidate = None
    sugar_reason = "No candidate deoxyribose sugar ring found."
    for ring in atom_rings:
        if len(ring) != 5:
            continue  # only consider 5-membered rings
        
        # Count ring oxygen atoms and carbon atoms.
        o_count = 0
        c_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "O":
                o_count += 1
            elif atom.GetSymbol() == "C":
                c_count += 1
            else:
                # If there are other atoms in the ring, skip this ring.
                o_count = -1
                break
        if o_count != 1 or c_count != 4:
            continue  # not a furanose-like ring
        
        # Next, count hydroxyl substituents on ring carbon atoms.
        # We consider a substituent to be a hydroxyl if it is a non-ring oxygen atom with exactly one neighbor (the ring carbon).
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "O":
                    # Check if this oxygen is a hydroxyl: degree==2 in mol (one from the carbon, one hydrogen)
                    # (Because we added explicit hydrogens, a typical –OH will have 2 neighbors.)
                    if nbr.GetDegree() == 2:
                        # further check that at least one neighbor is hydrogen
                        hasH = any(n.GetSymbol() == "H" for n in nbr.GetNeighbors())
                        if hasH:
                            hydroxyl_count += 1
        # In ribose, typically there are 3 OH groups attached to ring carbons;
        # in deoxyribose one OH (at the 2'-position) is missing, so we expect 2.
        if hydroxyl_count == 2:
            sugar_candidate = set(ring)
            sugar_reason = "Found candidate deoxyribose sugar ring (5-membered, 1 O + 4 C, 2 hydroxyl substituents)."
            break
    if sugar_candidate is None:
        return False, sugar_reason

    # STEP 2: Find candidate pyrimidine base rings.
    # We define a pyrimidine base as an aromatic, six-membered ring with exactly two nitrogen atoms.
    base_candidate = None
    base_reason = "No candidate pyrimidine base found."
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        # Check if all atoms are aromatic.
        all_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not all_aromatic:
            continue
        # Count nitrogen atoms in the ring.
        n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "N")
        if n_count != 2:
            continue
        base_candidate = set(ring)
        base_reason = "Found candidate pyrimidine base (6-membered aromatic ring with 2 nitrogen atoms)."
        break
    if base_candidate is None:
        return False, base_reason

    # STEP 3: Verify that the sugar and the base are connected via a glycosidic bond.
    # Look for any bond that connects an atom in the sugar ring to an atom in the base ring.
    glyco_bond_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (a1 in sugar_candidate and a2 in base_candidate) or (a2 in sugar_candidate and a1 in base_candidate):
            glyco_bond_found = True
            break
    if not glyco_bond_found:
        return False, "Sugar and base not connected via a glycosidic bond."

    # If all conditions are met, classify as pyrimidine deoxyribonucleoside.
    return True, "Molecule contains a deoxyribose sugar (5-membered ring with 2 hydroxyl substituents) connected to a pyrimidine base (6-membered aromatic ring with 2 nitrogens)."

# (Optional) A main section for quick testing.
if __name__ == '__main__':
    # Example: thymidine SMILES as given in the problem statement.
    examples = [
        "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O",  # thymidine
        "OC[C@H]1O[C@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O"          # 2'-deoxyuridine
    ]
    for smi in examples:
        res, reason = is_pyrimidine_deoxyribonucleoside(smi)
        print(f"SMILES: {smi}\nResult: {res}, Reason: {reason}\n")