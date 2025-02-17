"""
Classifies: CHEBI:28494 cardiolipin
"""
"""
Classifies: cardiolipin
A cardiolipin is defined as a phosphatidylglycerol composed of two molecules of phosphatidic acid
covalently linked to a molecule of glycerol. Thus, a typical cardiolipin should have:
  - Exactly 2 phosphorus atoms (each from one phosphatidate moiety),
  - At least 4 acyl ester groups (OC(=O)) indicating the 4 fatty acyl chains,
  - A glycerol linker connecting two phosphate groups. In many cardiolipins the phosphate–oxygen that
    connects to the glycerol appears in pairs that are separated by a path of 4 bonds (O–C–C–C–O)
    where the three middle atoms are carbons.
    
This improved heuristic is not perfect but attempts to reduce false positives by ensuring that the 
linkage between the phosphates is consistent with a glycerol backbone.
"""

from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on heuristics applied to its SMILES string.
        
    The heuristics are:
      1. There must be exactly 2 phosphorus atoms in the molecule.
      2. There should be at least 4 acyl ester bonds (pattern "OC(=O)") for the fatty acyl chains.
      3. There must be evidence for a glycerol linker connecting the two phosphatidate groups.
         This is checked by finding at least one pair of oxygen atoms (each directly attached to a phosphorus,
         but from different phosphate groups) that are connected by a shortest path of exactly 4 bonds
         where the 3 intermediate atoms are all carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a cardiolipin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for exactly 2 phosphorus atoms.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 2:
        return False, f"Expected 2 phosphorus atoms, found {len(phosphorus_atoms)}"
    
    # 2. Check for at least 4 acyl ester groups.
    # Use a SMARTS pattern for "OC(=O)"
    ester_pattern = Chem.MolFromSmarts("O-C(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Expected at least 4 acyl ester groups, found {len(ester_matches)}"
    
    # 3. Look for evidence of a glycerol linker connecting the two phosphatidate groups.
    # We will search for oxygen atoms that are directly attached to phosphorus.
    # For each oxygen that is directly attached to a phosphorus, record its index and the phosphorus index.
    op_oxygens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 15:
                    # record this oxygen along with the phosphorus atom index it is attached to
                    op_oxygens.append((atom.GetIdx(), nbr.GetIdx()))
                    break  # only need to record once per oxygen

    # We need to have at least one pair of such oxygens coming from different phosphorus atoms that are linked
    # via a glycerol backbone. Our heuristic: the shortest path from one oxygen to the other should span 5 atoms:
    # O -- C -- C -- C -- O, i.e. 4 bonds, and the three middle atoms must be carbons.
    found_linker = False
    n = len(op_oxygens)
    for i in range(n):
        for j in range(i+1, n):
            o1, p1 = op_oxygens[i]
            o2, p2 = op_oxygens[j]
            # Only consider oxygens originating from different phosphorus atoms.
            if p1 == p2:
                continue
            path = Chem.GetShortestPath(mol, o1, o2)
            # We expect a path that covers O-C-C-C-O (5 atoms, 4 bonds).
            if len(path) == 5:
                # Check that the 3 intermediate atoms (path[1], path[2], path[3]) are all carbons.
                if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in path[1:-1]):
                    found_linker = True
                    break
        if found_linker:
            break

    if not found_linker:
        return False, "Glycerol linker not detected in the headgroup (no oxygen pair connected via O–C–C–C–O path found)"
    
    # All heuristics pass: classify as cardiolipin.
    return True, ("Contains exactly 2 phosphorus atoms, at least 4 acyl ester groups, "
                  "and a glycerol linker (evidenced by an O–C–C–C–O pattern connecting "
                  "phosphorus-bound oxygens) consistent with cardiolipin.")

# Example usage (for testing purposes)
if __name__ == "__main__":
    # An example SMILES string for cardiolipin (one of the provided examples)
    example_smiles = ("P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)"
                      "(OC[C@H](O)COP(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)"
                      "CCCCCCC/C=C\\C/C=C\\CCCCC)(O)=O)(O)=O")
    result, reason = is_cardiolipin(example_smiles)
    print("Is cardiolipin?", result)
    print("Reason:", reason)