"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
#!/usr/bin/env python
"""
Classifies: pyrimidine deoxyribonucleoside

A pyrimidine deoxyribonucleoside is defined as a nucleoside that contains a deoxyribose sugar 
(i.e., lacking a 2'-OH) linked via a glycosidic bond to a pyrimidine base. In our heuristic, we:
  1. Reject molecules that contain any phosphorus (P) atoms (to avoid nucleotides and larger P–containing compounds).
  2. Look through rings to find a candidate deoxyribose sugar. We require that the sugar ring:
       • is five–membered (a furanose) made up of exactly four carbon atoms and one oxygen.
       • shows one –OH substituent attached to a ring carbon (since in deoxyribose the exocyclic CH2OH is not in the ring,
         and the 2'-OH is missing).
  3. Look for a candidate pyrimidine base that is a six–membered aromatic ring containing exactly two nitrogen atoms.
       Additionally, we require that the base carries at least one carbonyl (C=O), a common feature in pyrimidine bases.
  4. Finally, require that there is at least one bond connecting an atom in the sugar candidate to an atom in the base candidate.
  
If any check fails, a suitable explaining message is returned.
"""

from rdkit import Chem

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a pyrimidine deoxyribonucleoside, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Exclude molecules containing phosphorus (likely nucleotides or larger derivatives).
    if any(atom.GetSymbol() == "P" for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus atoms and is likely not a nucleoside."
    
    # To help detect –OH groups, add explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    sugar_candidate = None
    sugar_reason = "No candidate deoxyribose sugar ring found."
    # STEP 1: Search for a candidate sugar ring.
    for ring in atom_rings:
        # We expect a furanose ring: 5 members (4 carbons, 1 oxygen).
        if len(ring) != 5:
            continue
        # Count number of carbons and oxygens in the ring.
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
        
        # Count hydroxyl substituents on the ring carbons.
        # (We only consider neighbors not in the ring that are oxygen atoms and have at least one hydrogen attached.)
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != "C":
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "O":
                    # Check if this oxygen has any H neighbor.
                    if any(n.GetSymbol() == "H" for n in nbr.GetNeighbors()):
                        hydroxyl_count += 1
        # For a deoxyribose, the only OH on the sugar ring should be that on the 3'-carbon.
        if hydroxyl_count == 1:
            sugar_candidate = set(ring)
            sugar_reason = "Found candidate deoxyribose sugar ring (5-membered, 4 C + 1 O, 1 hydroxyl on ring)."
            break
    if sugar_candidate is None:
        return False, sugar_reason
    
    # STEP 2: Find candidate pyrimidine base rings.
    base_candidate = None
    base_reason = "No candidate pyrimidine base found."
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        # Check that every atom in the ring is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count number of nitrogen atoms.
        n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "N")
        if n_count != 2:
            continue
        
        # Additionally, verify that at least one ring atom (or its immediate neighbor outside the ring)
        # carries a carbonyl (C=O). Look over all bonds around ring atoms.
        carbonyl_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider carbon atoms inside the ring.
            if atom.GetSymbol() != "C":
                continue
            for bond in atom.GetBonds():
                # If bond is double and to oxygen and the oxygen is not in the ring, consider it a carbonyl.
                if bond.GetBondTypeAsDouble() == 2.0:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetSymbol() == "O" and nbr.GetIdx() not in ring:
                        carbonyl_found = True
                        break
            if carbonyl_found:
                break
        if not carbonyl_found:
            continue
        
        base_candidate = set(ring)
        base_reason = "Found candidate pyrimidine base (6-membered aromatic ring with 2 nitrogens and at least one carbonyl)."
        break
    if base_candidate is None:
        return False, base_reason
    
    # STEP 3: Check for a glycosidic bond between the sugar candidate and the base candidate.
    glyco_bond_found = False
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if (a1 in sugar_candidate and a2 in base_candidate) or (a2 in sugar_candidate and a1 in base_candidate):
            glyco_bond_found = True
            break
    if not glyco_bond_found:
        return False, "Sugar and base are not connected via a glycosidic bond."
    
    return True, "Molecule contains a deoxyribose sugar (5-membered, 4 C + 1 O, 1 hydroxyl on ring) connected to a pyrimidine base (6-membered aromatic ring with 2 N and carbonyl)."

# (Optional) Quick testing when run as a script.
if __name__ == '__main__':
    # Examples taken from the problem statement.
    examples = [
        "Cc1cn([C@H]2C[C@H](O)[C@@H](CO)O2)c(=O)[nH]c1=O",  # thymidine
        "OC[C@H]1O[C@H](C[C@@H]1O)n1ccc(=O)[nH]c1=O"         # 2'-deoxyuridine
    ]
    for smi in examples:
        res, reason = is_pyrimidine_deoxyribonucleoside(smi)
        print(f"SMILES: {smi}\nResult: {res}, Reason: {reason}\n")