"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar – any monosaccharide (or nucleoside) containing an alcoholic –OH 
esterified with phosphoric acid.
Definition: The molecule must contain at least one candidate sugar ring (a 5-membered furanose
or 6-membered pyranose with one ring oxygen, sufficient carbons, and at least two exocyclic –OH groups)
or an open-chain sugar fragment (a contiguous polyol pattern), and at least one phosphate group 
(a phosphorus atom with at least one double‐bonded oxygen) that is linked via an alcoholic –O to a sugar–derived carbon.
"""

from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar by first searching for a candidate sugar ring
    (and confirming that it bears at least two exocyclic hydroxyl groups) or, if no ring is found,
    by looking for an open-chain sugar fragment. Then the molecule is confirmed only if at least one phosphate 
    group (P atom with one double-bonded O) is attached via a single-bonded oxygen to a carbon associated with the sugar.

    Args:
        smiles (str): SMILES string for the molecule.

    Returns:
        bool: True if classified as phospho sugar; False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Identify candidate sugar rings based on ring size, composition and exocyclic –OH count.
    # A candidate sugar ring is either:
    #    - 5-membered: 1 oxygen and 4 carbons, OR
    #    - 6-membered: 1 oxygen and 5 carbons.
    # Also require that at least 2 exocyclic –OH groups are attached to ring carbons.
    candidate_sugar_rings = []   # each ring represented as a set of atom indices
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
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                count_O += 1
            elif atomic_num == 6:
                count_C += 1
            elif atomic_num == 15:
                contains_P = True
        # Skip rings that contain P or do not have exactly 1 oxygen and the right number of carbons
        if contains_P:
            continue
        if (size == 5 and (count_O != 1 or count_C != 4)) or (size == 6 and (count_O != 1 or count_C != 5)):
            continue
        
        # Now check for exocyclic hydroxyls attached to ring carbons.
        hx_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look for oxygen neighbors not in the ring that appear as hydroxyl groups
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue  # skip in-ring neighbors
                if nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() > 0:
                    hx_count += 1
        if hx_count >= 2:
            candidate_sugar_rings.append(set(ring))
    
    # STEP 2: Identify phosphate groups.
    # A phosphate group is a phosphorus atom (atomic #15) that has at least one double-bonded oxygen.
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
    
    # STEP 3: For each phosphate, check if any of the oxygen substituents (attached via a single bond)
    # is connected (through a carbon) to a sugar ring or its immediate substituent.
    for p_atom in phosphate_atoms:
        for o_atom in p_atom.GetNeighbors():
            if o_atom.GetAtomicNum() != 8:
                continue  # must be oxygen
            bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), o_atom.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            # Now inspect the neighbors of this phosphate-bound oxygen (except the phosphorus itself)
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetIdx() == p_atom.GetIdx():
                    continue
                if neighbor.GetAtomicNum() == 6:  # carbon neighbor
                    cid = neighbor.GetIdx()
                    # Test if the carbon is directly in a candidate sugar ring …
                    if any(cid in ring for ring in candidate_sugar_rings):
                        return True, ("Found candidate sugar ring with sufficient hydroxyls and a phosphate group "
                                      "attached via an alcoholic O to a sugar-derived carbon")
                    # … or if one of its neighbors is in a sugar ring.
                    for nb in neighbor.GetNeighbors():
                        if any(nb.GetIdx() in ring for ring in candidate_sugar_rings):
                            return True, ("Found a sugar ring substituent (carbon adjacent to a candidate ring) with phosphate attachment")
    
    # STEP 4: Fallback for open-chain sugar fragments.
    # Some sugars in non-cyclic form display a contiguous polyol pattern.
    # We use a SMARTS that looks for three consecutive carbons each bearing an OH.
    sugar_chain_smarts = Chem.MolFromSmarts("C(O)C(O)C(O)")
    if sugar_chain_smarts is not None:
        chain_matches = mol.GetSubstructMatches(sugar_chain_smarts)
        if chain_matches:
            # For each phosphate group, check if its alcoholic oxygen attaches through a carbon that is part of a sugar chain.
            for p_atom in phosphate_atoms:
                for o_atom in p_atom.GetNeighbors():
                    if o_atom.GetAtomicNum() != 8:
                        continue
                    bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), o_atom.GetIdx())
                    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                        continue
                    for neighbor in o_atom.GetNeighbors():
                        if neighbor.GetIdx() == p_atom.GetIdx():
                            continue
                        if neighbor.GetAtomicNum() == 6:
                            for match in chain_matches:
                                if neighbor.GetIdx() in match:
                                    return True, ("Found an open-chain sugar fragment (polyol pattern) with phosphate attached "
                                                   "via an alcoholic O to a sugar-derived carbon")
    
    return False, "No alcoholic –OH on a sugar ring or open-chain sugar fragment found that is esterified to phosphate"

# When run as a script, several test examples are provided.
if __name__ == "__main__":
    test_examples = [
        # True positives:
        ("P(OCC1O[C@H](O)C(O)C(O)C1O)(O)(O)=O", "alpha-D-Hexose 6-phosphate"),
        ("Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N", "2'-deoxy-5-methyl-5'-cytidylic acid"),
        ("[C@@H]1O[C@H](OP(O)(O)=O)[C@H](O)[C@@H]1O", "alpha-L-fucose 1-phosphate"),
        # A false positive example from the previous attempt:
        ("C1[C@H]2[C@@H]([C@@H]([C@H](O2)N3C4=C(C(=NC=N4)N)N=C3SC5=CC=C(C=C5)Cl)O)OP(=O)(O1)O", 
         "False positive example (should not be classified)"),
        # An example that was missed previously (open-chain sugar form):
        ("P(OC[C@@H](O)[C@@H](O)[C@H](O)C(=O)CO)([O-])([O-])=O", "D-Fructose 6-Phosphate-Disodium Salt"),
    ]
    
    for smi, name in test_examples:
        result, explanation = is_phospho_sugar(smi)
        print(f"Name: {name}\nSMILES: {smi}\nClassified as phospho sugar? {result}\nReason: {explanation}\n")