"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: cardiolipin
Cardiolipin is defined as a phosphatidylglycerol composed of two molecules of phosphatidic acid 
covalently linked to a molecule of glycerol. Heuristics used:
  1. The molecule must be a single connected structure (no counterions or additional fragments).
  2. It must have an overall formal charge of 0 (to avoid misclassification of ionized forms).
  3. It must contain exactly 2 phosphorus atoms.
  4. It must contain at least 4 acyl ester groups, represented by an "O-C(=O)" substructure.
  5. It must show evidence of a glycerol linker connecting the two phosphatidate groups,
     detected by a shortest path of exactly 5 atoms (O–C–C–C–O) between phosphorus-bound 
     oxygen atoms (excluding those that are part of an acyl ester).
"""

from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string using several heuristics.
    
    Steps:
      1. Ensure the SMILES represents a single connected, neutral molecule.
      2. Check that the overall formal charge of the molecule is 0.
      3. Verify that the molecule contains exactly 2 phosphorus atoms.
      4. Find at least 4 acyl ester groups (substructure "O-C(=O)").
      5. Identify oxygen atoms directly bonded to a phosphorus that are NOT part of an acyl ester.
         Then, check if at least one pair of these oxygens (coming from different phosphorus atoms)
         is connected by a path of exactly 5 atoms (O–C–C–C–O) where the middle three atoms are carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as cardiolipin, False otherwise.
        str: Explanation for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check connectivity: The molecule should be a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=False)
    if len(frags) != 1:
        return False, f"Expected 1 connected component, found {len(frags)} fragments (possible counterions present)"
        
    # 2. Check overall formal charge is 0.
    if Chem.GetFormalCharge(mol) != 0:
        return False, f"Overall formal charge ({Chem.GetFormalCharge(mol)}) is non-zero; ionized species are not accepted"
    
    # 3. Count exactly 2 phosphorus atoms.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 2:
        return False, f"Expected 2 phosphorus atoms, found {len(phosphorus_atoms)}"
    
    # 4. Check for at least 4 acyl ester groups: using the SMARTS "O-C(=O)".
    ester_pattern = Chem.MolFromSmarts("O-C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Expected at least 4 acyl ester groups, found {len(ester_matches)}"
    
    # 5. Identify oxygen atoms directly bonded to phosphorus but NOT part of an acyl ester.
    #    We will collect indices of those oxygens along with the index of the phosphorus to which they are attached.
    op_oxygens = []  # list of tuples (oxygen_atom_idx, phosphorus_atom_idx)
    for oxygen in mol.GetAtoms():
        if oxygen.GetAtomicNum() != 8:
            continue
        # Check if oxygen is bonded to a phosphorus.
        p_idx = None
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() == 15:
                p_idx = nbr.GetIdx()
                break
        if p_idx is None:
            continue
        
        # Now check if this oxygen is part of an acyl ester.
        is_ester = False
        for nbr in oxygen.GetNeighbors():
            if nbr.GetAtomicNum() == 15:
                continue  # skip the phosphorus neighbor
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                # Look through bonds of the carbon to check for a carbonyl double bond (C=O).
                for bond in nbr.GetBonds():
                    # bond.GetBondTypeAsDouble() returns 2.0 for a double bond.
                    if bond.GetBondTypeAsDouble() == 2:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8 and other.GetIdx() != oxygen.GetIdx():
                            is_ester = True
                            break
                if is_ester:
                    break
        if not is_ester:
            op_oxygens.append((oxygen.GetIdx(), p_idx))
    
    if len(op_oxygens) < 2:
        return False, "Insufficient phosphorus-bound non-ester oxygens to detect a glycerol linker"
    
    # 6. Look for at least one pair of these oxygen atoms (from different P atoms)
    #    connected by a path of exactly 5 atoms (indicative of a glycerol linker: O–C–C–C–O).
    found_linker = False
    n = len(op_oxygens)
    for i in range(n):
        for j in range(i+1, n):
            o1_idx, p1_idx = op_oxygens[i]
            o2_idx, p2_idx = op_oxygens[j]
            if p1_idx == p2_idx:
                continue  # must come from different phosphorus atoms
            path = Chem.GetShortestPath(mol, o1_idx, o2_idx)
            # We require a path consisting of exactly 5 atoms (4 bonds): O–C–C–C–O.
            if len(path) == 5:
                # Check that the three intermediate atoms are all carbons.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path[1:-1]):
                    found_linker = True
                    break
        if found_linker:
            break
    
    if not found_linker:
        return False, "Glycerol linker not detected: no O–C–C–C–O path found connecting phosphorus-bound non-ester oxygens"
    
    return True, ("Contains exactly 2 phosphorus atoms, at least 4 acyl ester groups, and a glycerol linker (via an O–C–C–C–O path between "
                  "phosphorus-bound non-ester oxygens) consistent with cardiolipin.")

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        # A correct neutral cardiolipin example.
        "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OC[C@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O)(O)=O",
        # An ionized (false positive) example should fail due to a non-zero formal charge.
        "CCCCCCCC\\C=C/CCCCCCCC(=O)O[C@H](COP([O-])(=O)OCC(O)COP([O-])(=O)OC[C@@H](COC(=O)CCCCCCC\\C=C/CCCCCCCC)OC(=O)CCCCCCC\\C=C/CCCCCCCC)OC(=O)CCCCCCC\\C=C/CCCCCCCC"
    ]
    
    for smi in test_smiles:
        result, reason = is_cardiolipin(smi)
        print("SMILES:", smi)
        print("Is cardiolipin?", result)
        print("Reason:", reason)
        print("-" * 80)