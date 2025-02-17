"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: cardiolipin
Cardiolipin is defined as a phosphatidylglycerol composed of two molecules of phosphatidic acid 
covalently linked to a molecule of glycerol. Heuristics used:
  1. The molecule must be a single connected structure (no counterions or additional fragments).
  2. It must contain exactly 2 phosphorus atoms.
  3. It must contain at least 4 acyl ester groups, represented by the "O–C(=O)" substructure.
  4. It must show evidence of a glycerol linker connecting the two phosphatidate groups. This is 
     detected by finding at least one pair of phosphorus-bound oxygen atoms (excluding those part 
     of an acyl ester) that are linked by a shortest path of exactly 5 atoms (O–C–C–C–O) in which 
     the three inner atoms are all carbons.
"""

from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string using several heuristics.
    
    Steps:
      1. Ensure the SMILES represents a single connected molecule.
      2. Count exactly 2 phosphorus atoms (atomic number 15) in the molecule.
      3. Count the number of acyl ester groups via the SMARTS pattern "O-C(=O)" and ensure at least 4.
      4. Identify oxygen atoms that are directly attached to phosphorus but NOT part of an acyl ester.
         Then, check if there exists at least one pair (coming from different phosphorus atoms) that 
         are connected by a path exactly 5 atoms long (O–C–C–C–O) where the 3 intervening atoms are carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as cardiolipin, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check molecule connectivity: it should be a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=False)
    if len(frags) != 1:
        return False, f"Expected 1 connected component, found {len(frags)} fragments (possible counterions present)"
    
    # 2. Verify exactly 2 phosphorus atoms.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 2:
        return False, f"Expected 2 phosphorus atoms, found {len(phosphorus_atoms)}"
    
    # 3. Check for at least 4 acyl ester groups.
    ester_pattern = Chem.MolFromSmarts("O-C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Expected at least 4 acyl ester groups, found {len(ester_matches)}"
    
    # 4. Identify oxygen atoms directly bound to phosphorus 
    #     but not part of an acyl ester. For each oxygen, we check if it is bonded to a phosphorus,
    #     and then ignore it if one of its other bonds (to a carbon) is part of a carbonyl (C=O).
    op_oxygens = []  # list of tuples (oxygen_atom_idx, phosphorus_atom_idx)
    for oxygen in mol.GetAtoms():
        if oxygen.GetAtomicNum() != 8:
            continue
        # Check if oxygen is bonded to a phosphorus atom.
        p_at_idx = None
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() == 15:
                p_at_idx = nbr.GetIdx()
                break
        if p_at_idx is None:
            continue  # oxygen not bound to P
        
        # Check if this oxygen is part of an acyl ester.
        # That is, if it is also bound (besides the P) to a carbon that has a double bond to an oxygen.
        is_ester = False
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() == 15:
                continue  # skip the phosphorus neighbor
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                # check for a double bond from this carbon to an oxygen (i.e. a carbonyl)
                for bond in nbr.GetBonds():
                    # if the bond is a double bond and the other atom is oxygen (and not our original oxygen)
                    if bond.GetBondTypeAsDouble() == 2:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8 and other.GetIdx() != oxygen.GetIdx():
                            is_ester = True
                            break
                if is_ester:
                    break
        if not is_ester:
            op_oxygens.append((oxygen.GetIdx(), p_at_idx))
    
    if len(op_oxygens) < 2:
        return False, "Insufficient phosphorus-bound non-ester oxygens to detect a glycerol linker"
    
    # 5. Look for at least one pair of these oxygen atoms (attached to different phosphorus atoms)
    #    that are connected by a path: O – C – C – C – O (5 atoms, 4 bonds) with the three inner atoms all carbons.
    found_linker = False
    n = len(op_oxygens)
    for i in range(n):
        for j in range(i+1, n):
            o1, p1 = op_oxygens[i]
            o2, p2 = op_oxygens[j]
            if p1 == p2:
                continue  # must come from different phosphorus atoms
            path = Chem.GetShortestPath(mol, o1, o2)
            # We expect a path that covers 5 atoms (4 bonds): O-C-C-C-O.
            if len(path) == 5:
                # Check that the middle three atoms are all carbons.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path[1:-1]):
                    found_linker = True
                    break
        if found_linker:
            break
    
    if not found_linker:
        return False, "Glycerol linker not detected (no O–C–C–C–O path found connecting phosphorus-bound non-ester oxygens)"
    
    return True, ("Contains exactly 2 phosphorus atoms, at least 4 acyl ester groups, and a glycerol linker (evidenced by an O–C–C–C–O pattern connecting "
                  "phosphorus-bound non-ester oxygens) consistent with cardiolipin.")

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        # A correct cardiolipin example from the provided list.
        "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OC[C@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O)(O)=O",
        # A false positive (salt) example – note the dot (".") indicates additional fragments.
        "[NH4+].[NH4+].CCCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@H](O)COP([O-])(=O)OC[C@@H](COC(=O)CCCCCCCCCCCN)OC(=O)CCCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC"
    ]
    
    for smi in test_smiles:
        result, reason = is_cardiolipin(smi)
        print("SMILES:", smi)
        print("Is cardiolipin?", result)
        print("Reason:", reason)
        print("-" * 80)