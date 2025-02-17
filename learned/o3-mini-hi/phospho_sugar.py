"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar – any monosaccharide (or nucleoside) containing an alcoholic –OH 
esterified with phosphoric acid.

Definition: The molecule must contain at least one candidate sugar ring (a 5-membered ring 
with exactly 4 carbons and 1 oxygen, or a 6-membered ring with 5 carbons and 1 oxygen) 
and at least one phosphate group (a phosphorus atom having at least one double-bonded oxygen) 
that is linked via a single-bonded oxygen to a carbon directly in, or immediately adjacent 
to, such a ring.
"""

from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    
    The detection process is as follows:
      1. Identify candidate sugar rings from ring information. A candidate sugar ring is:
           - Either 5-membered (furanose: 4 carbons and 1 oxygen) or 6-membered (pyranose:
             5 carbons and 1 oxygen).
           - The ring should not contain phosphorus.
      2. Identify phosphate groups:
           - Look for phosphorus atoms (atomic number 15) with at least one double‐bonded oxygen.
      3. For each phosphate group, look at all oxygen substituents attached via a single bond.
         Then for each such oxygen, check its non–phosphorus neighbor carbons. If one of these
         carbons is either part of a candidate sugar ring or is directly attached (neighbor) to an 
         atom in a candidate sugar ring, the molecule is classified as a phospho sugar.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phospho sugar, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Identify candidate sugar rings based on ring size and atomic composition.
    candidate_sugar_rings = []  # list of sets, each set is the atom indices in a qualifying ring
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        size = len(ring)
        if size not in (5, 6):
            continue
        count_O = 0
        count_C = 0
        contains_P = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            num = atom.GetAtomicNum()
            if num == 8:
                count_O += 1
            elif num == 6:
                count_C += 1
            elif num == 15:
                contains_P = True
            # We ignore other atoms for simplicity.
        # For a furanose: expect size==5 with 1 oxygen and 4 carbons.
        # For a pyranose: expect size==6 with 1 oxygen and 5 carbons.
        if contains_P:
            continue
        if size == 5 and count_O == 1 and count_C == 4:
            candidate_sugar_rings.append(set(ring))
        elif size == 6 and count_O == 1 and count_C == 5:
            candidate_sugar_rings.append(set(ring))
            
    if not candidate_sugar_rings:
        return False, "No candidate sugar ring (5- or 6-membered with expected C/O count) found"
    
    # STEP 2: Identify phosphate groups: phosphorus with at least one double-bonded oxygen.
    phosphate_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        has_double_O = False
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    has_double_O = True
                    break
        if has_double_O:
            phosphate_atoms.append(atom)
    
    if not phosphate_atoms:
        return False, "No phosphate group (P with at least one double-bonded O) found"
    
    # STEP 3: For each phosphate group, check its oxygen substituents (via a single bond)
    # to see if they attach to a carbon that is part of, or directly adjacent to, a candidate sugar ring.
    for p_atom in phosphate_atoms:
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetAtomicNum() != 8:
                continue  # Must be oxygen
            bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), o_atom.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # For the oxygen attached via a single bond, inspect its other neighbors.
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetIdx() == p_atom.GetIdx():
                    continue
                # Look for a carbon neighbor
                if neighbor.GetAtomicNum() == 6:
                    carbon_idx = neighbor.GetIdx()
                    # Check if this carbon is directly in any candidate sugar ring
                    in_ring = any(carbon_idx in ring for ring in candidate_sugar_rings)
                    # OR check if any neighbor of this carbon is in a candidate sugar ring
                    adjacent_to_ring = any(any(nb.GetIdx() in ring for ring in candidate_sugar_rings)
                                            for nb in neighbor.GetNeighbors())
                    if in_ring or adjacent_to_ring:
                        return True, ("Found candidate sugar ring (or its immediate substituent) "
                                      "and a phosphate group attached via an alcoholic O to a sugar-derived C")
    
    return False, "No alcoholic –OH on a candidate sugar ring (or immediate substituent) found that is esterified to a phosphate"

# For testing the function (executed when the module is run as a script)
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("P(OCC1O[C@H](O)C(O)C(O)C1O)(O)(O)=O", "alpha-D-Hexose 6-phosphate"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N", "2'-deoxy-5-methyl-5'-cytidylic acid"),
        ("[C@@H]1O[C@H](OP(O)(O)=O)[C@H](O)[C@@H]1O", "alpha-L-fucose 1-phosphate"),
        ("C1[C@H]2[C@@H]([C@@H]([C@H](O2)N3C4=C(C(=NC=N4)N)N=C3SC5=CC=C(C=C5)Cl)O)OP(=O)(O1)O", 
         "False positive example (should not be classified)"),
    ]
    
    for smi, name in test_examples:
        classified, reason = is_phospho_sugar(smi)
        print(f"Name: {name}\nSMILES: {smi}\nClassified as phospho sugar? {classified}\nReason: {reason}\n")