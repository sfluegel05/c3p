"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar – any monosaccharide (or nucleoside) containing an alcoholic –OH 
esterified with phosphoric acid.
Definition: A molecule must contain at least one small sugar ring (5‐ or 6‐membered ring 
with exactly one ring oxygen) and at least one phosphate group in which one oxygen (bonded 
via a single bond to phosphorus) is connected (through a carbon) to that ring.
"""

from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    
    The detection is implemented as follows:
    1. The algorithm finds candidate sugar rings by extracting all rings (from mol.GetRingInfo)
       that are of size 5 or 6 and contain exactly one oxygen atom. This typical signature is used
       for furanose or pyranose rings.
    2. It then looks for phosphate groups – defined as any phosphorus atom (atomic number 15)
       that has at least one double-bonded oxygen.
    3. For every phosphate group found, it examines each neighbor oxygen (bonded via a single 
       bond) and then checks whether that oxygen is attached (via a carbon atom) to any candidate
       sugar ring.
       
    Args:
        smiles (str): SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phospho sugar, False otherwise.
        str: A textual explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify candidate sugar rings:
    # We use the ring info from RDKit and filter for rings of size 5 or 6 that have exactly one oxygen.
    sugar_rings = []
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) in [5, 6]:
            o_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if o_count == 1:
                sugar_rings.append(set(ring))
    
    if not sugar_rings:
        return False, "No sugar ring (size 5 or 6 with one oxygen) found"
    
    # Identify phosphate groups:
    # A phosphate group is represented here by a phosphorus (atomic number 15) that has at least 
    # one double-bonded oxygen neighbor.
    phosphate_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus
            has_dbl_O = False
            for nb in atom.GetNeighbors():
                if nb.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                        has_dbl_O = True
                        break
            if has_dbl_O:
                phosphate_atoms.append(atom)
    
    if not phosphate_atoms:
        return False, "No phosphate group (P with at least one =O bond) found"
    
    # Now, for each phosphate group, check its oxygen substituents.
    # We want to see if any oxygen (attached by a single bond to phosphorus) bridges to a carbon that belongs to a sugar ring.
    for p_atom in phosphate_atoms:
        for nb in p_atom.GetNeighbors():
            if nb.GetAtomicNum() != 8:
                continue  # skip non-oxygen neighbors
            bond = mol.GetBondBetweenAtoms(p_atom.GetIdx(), nb.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue  # we only care for P–O single bonds
            # For each oxygen bonded to P, check its neighbors (besides the phosphorus)
            for oxy_nb in nb.GetNeighbors():
                if oxy_nb.GetIdx() == p_atom.GetIdx():
                    continue
                # Look for a carbon neighbor that is part of one of our candidate sugar rings.
                if oxy_nb.GetAtomicNum() == 6:
                    for ring in sugar_rings:
                        if oxy_nb.GetIdx() in ring:
                            return True, ("Found sugar ring (size 5 or 6 with one oxygen) "
                                          "with an alcoholic –OH (via a bridging oxygen) esterified to a phosphoric acid")
    return False, "No alcoholic –OH on a candidate sugar ring found that is esterified to phosphoric acid"


# For testing the function (execute when running as script):
if __name__ == "__main__":
    # Test examples: (feel free to add more tests)
    test_smiles = [
        "P(OCC1O[C@H](O)C(O)C(O)C1O)(O)(O)=O",  # alpha-D-Hexose 6-phosphate ~ expected True
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N",  # 2'-deoxy-5-methyl-5'-cytidylic acid ~ expected True
        # A false positive example: a molecule that has a sugar ring and a phosphate group,
        # but the phosphate is not connected to the sugar ring (this is only an example string)
        "C1[C@H]2[C@@H]([C@@H]([C@H](O2)N3C4=C(C(=NC=N4)N)N=C3SC5=CC=C(C=C5)Cl)O)OP(=O)(O1)O"
    ]
    
    for smi in test_smiles:
        classification, reason = is_phospho_sugar(smi)
        print(f"SMILES: {smi}\nClassified as phospho sugar? {classification}\nReason: {reason}\n")