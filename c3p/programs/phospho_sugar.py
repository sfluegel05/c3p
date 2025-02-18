"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar – any monosaccharide (or nucleoside) containing an alcoholic –OH 
esterified with phosphoric acid.

Definition: The molecule must contain at least one sugar ring (furanose or pyranose) 
with the expected atomic makeup – a 5-membered ring with 4 carbons and 1 oxygen (furanose) 
or a 6-membered ring with 5 carbons and 1 oxygen (pyranose) – and at least one phosphate group 
(where a phosphorus atom is connected to at least one double-bonded oxygen). One of the non–P 
neighbors (via a single bond) of that phosphate must be an oxygen that is in turn bound to a carbon 
that is part of such a sugar ring.
"""

from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    
    The detection is implemented as follows:
      1. Identify all candidate sugar rings by examining the ring membership of atoms.
         For a candidate ring, we require:
           - Its size is 5 (furanose) or 6 (pyranose).
           - It has exactly one oxygen atom.
           - Its carbon count is 4 (if size is 5) or 5 (if size is 6).
           - None of the ring atoms is phosphorus.
    
      2. Identify phosphate groups as phosphorus atoms (atomic number 15) that have at 
         least one neighbor oxygen connected via a double bond.
    
      3. For each phosphate group, examine oxygen neighbors attached via a single bond.
         For each such oxygen, look at its other neighbor atoms. If one is a carbon that is 
         part of one of our candidate sugar rings then we classify the molecule as a phospho sugar.
    
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phospho sugar, False otherwise.
        str: A textual explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Gather candidate sugar rings from ring information.
    # We require: ring size 5 (4 C, 1 O) OR ring size 6 (5 C, 1 O) and no P in the ring.
    candidate_sugar_rings = []
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
            # (Other atoms are not expected in a simple monosaccharide ring)
        # For a furanose: 5-membered => 4 C's and 1 O; for a pyranose: 6-membered => 5 C's and 1 O.
        if contains_P:
            continue
        if size == 5 and count_O == 1 and count_C == 4:
            candidate_sugar_rings.append(set(ring))
        elif size == 6 and count_O == 1 and count_C == 5:
            candidate_sugar_rings.append(set(ring))
    
    if not candidate_sugar_rings:
        return False, "No candidate sugar ring (5- or 6-membered with expected C/O count) found"
    
    # Identify phosphate groups:
    # A phosphate group is indicated by a phosphorus atom (15) with at least one double-bonded oxygen.
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
    
    # For each phosphate group, check its oxygen substituents that are bound via a single bond.
    for p_atom in phosphate_atoms:
        for nb in p_atom.GetNeighbors():
            # Skip if neighbor is not oxygen.
            if nb.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), nb.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # For the oxygen attached to P via a single bond, check its other neighbor(s).
            for o_nb in nb.GetNeighbors():
                if o_nb.GetIdx() == p_atom.GetIdx():
                    continue
                # We need an oxygen attached to a carbon that belongs to a sugar ring.
                if o_nb.GetAtomicNum() == 6:
                    # Check if this carbon atom is in any candidate sugar ring.
                    for ring in candidate_sugar_rings:
                        if o_nb.GetIdx() in ring:
                            return True, ("Found candidate sugar ring (size 5 or 6 with expected C/O count) "
                                          "and a phosphate group attached via an alcoholic O leading to a ring C")
    return False, "No alcoholic –OH on a candidate sugar ring found that is esterified to a phosphate"

# For testing the function (executed when the module is run as a script)
if __name__ == "__main__":
    test_smiles = [
        "P(OCC1O[C@H](O)C(O)C(O)C1O)(O)(O)=O",  # alpha-D-Hexose 6-phosphate ~ expected True
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N",  # 2'-deoxy-5-methyl-5'-cytidylic acid ~ expected True
        "[C@@H]1O[C@H](OP(O)(O)=O)[C@H](O)[C@@H]1O",  # alpha-L-fucose 1-phosphate ~ expected True
        # A false positive example (phosphate attached to a fused ring system that does not meet our candidate sugar criteria)
        "C1[C@H]2[C@@H]([C@@H]([C@H](O2)N3C4=C(C(=NC=N4)N)N=C3SC5=CC=C(C=C5)Cl)O)OP(=O)(O1)O"
    ]
    
    for smi in test_smiles:
        classification, reason = is_phospho_sugar(smi)
        print(f"SMILES: {smi}\nClassified as phospho sugar? {classification}\nReason: {reason}\n")