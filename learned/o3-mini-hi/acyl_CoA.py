"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition:
  An acyl-CoA is a thioester that results from condensation of the thiol group of Coenzyme A
  with the carboxy group of a fatty acid. Its structure features an acyl thioester group
  (C(=O)S) directly attached to a pantetheine moiety (which in our case may be recognized by
  either of two SMARTS patterns: "SCCNC(=O)CCNC(=O)" or "NCCSC(=O)") that in turn is connected
  (via a phosphate bridge) to an adenine-containing nucleotide.
Improvement:
  Instead of only testing for the mere presence of three independent fragments, we:
    1. Look for an acyl thioester group "[CX3](=O)[SX2]".
    2. Look for the pantetheine fragment using two possible patterns.
    3. Look for an adenine moiety "n1cnc2c(N)ncnc12".
    4. Verify that at least one thioester sulfur is either part of a pantetheine match or
       has a neighboring atom belonging to a pantetheine match.
    5. Check that at least one pantetheine match is connected (by a short path, ≤25 bonds)
       to the adenine moiety with at least one phosphorus atom in that connecting path.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    The molecule must contain:
      1. An acyl thioester fragment: "[CX3](=O)[SX2]".
      2. A pantetheine fragment which we look for using two SMARTS:
         a) "SCCNC(=O)CCNC(=O)"
         b) "NCCSC(=O)"
      3. An adenine moiety: "n1cnc2c(N)ncnc12".
    In addition, we verify that:
      A. The sulfur in an acyl thioester is part of (or immediately adjacent to) at least one
         pantetheine match.
      B. At least one pantetheine match is connected to the adenine moiety by a short topological
         path (<=25 bonds) that includes at least one phosphorus (P) atom (i.e. representing a phosphate bridge).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecular features are consistent with an acyl-CoA.
        str: Explanation for the classification decision.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Locate the acyl thioester group.
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pat = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pat is None:
        return False, "Error creating thioester SMARTS"
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "Thioester group '[CX3](=O)[SX2]' not found"
    
    # 2. Locate the pantetheine fragment using two different SMARTS patterns and combine matches.
    pantetheine_patterns = []
    pantetheine_smarts1 = "SCCNC(=O)CCNC(=O)"
    pant_pat1 = Chem.MolFromSmarts(pantetheine_smarts1)
    if pant_pat1:
        pantetheine_patterns.append(pant_pat1)
    pantetheine_smarts2 = "NCCSC(=O)"
    pant_pat2 = Chem.MolFromSmarts(pantetheine_smarts2)
    if pant_pat2:
        pantetheine_patterns.append(pant_pat2)
    
    pantetheine_matches = []
    for pat in pantetheine_patterns:
        matches = mol.GetSubstructMatches(pat)
        if matches:
            pantetheine_matches.extend(matches)
    if not pantetheine_matches:
        return False, "Pantetheine fragment not detected using known patterns"
    
    # 3. Locate the adenine moiety.
    adenine_smarts = "n1cnc2c(N)ncnc12"
    adenine_pat = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pat is None:
        return False, "Error creating adenine SMARTS pattern"
    adenine_matches = mol.GetSubstructMatches(adenine_pat)
    if not adenine_matches:
        return False, "Adenine moiety of CoA not detected"
    
    # 4. Verify that a thioester sulfur is connected to the pantetheine fragment.
    # For each thioester match, index 1 corresponds to the sulfur.
    thioester_linked = False
    for tmatch in thioester_matches:
        thio_sulfur = tmatch[1]
        for pant_match in pantetheine_matches:
            # If the thioester sulfur is directly in the pantetheine match, that's good.
            if thio_sulfur in pant_match:
                thioester_linked = True
                break
            # Alternatively, if the sulfur has a neighbor that is in the pantetheine match.
            atom = mol.GetAtomWithIdx(thio_sulfur)
            nbr_idxs = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
            if any(n in pant_match for n in nbr_idxs):
                thioester_linked = True
                break
        if thioester_linked:
            break
    if not thioester_linked:
        return False, "Thioester group is not connected to the pantetheine fragment as expected"
    
    # 5. Check connectivity between a pantetheine match and an adenine moiety.
    max_path_length = 25
    found_connection = False
    for pant_match in pantetheine_matches:
        for adenine_match in adenine_matches:
            # For each combination of atoms from the two matches, get the shortest path.
            for pant_atom in pant_match:
                for aden_atom in adenine_match:
                    try:
                        path = Chem.GetShortestPath(mol, pant_atom, aden_atom)
                    except Exception:
                        continue
                    if path and len(path) <= max_path_length:
                        # Verify that there is at least one phosphorus (atomic number 15) along the path.
                        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 15 for idx in path):
                            found_connection = True
                            break
                if found_connection:
                    break
            if found_connection:
                break
        if found_connection:
            break
    if not found_connection:
        return False, "Pantetheine and adenine fragments are not connected via a proper phosphate linker (<=25 bonds with a phosphorus)"
    
    return True, "Contains acyl thioester attached to a pantetheine fragment and linked to an adenine moiety (acyl-CoA)"

# Example usage:
if __name__ == '__main__':
    # Testing with one of the provided acyl-CoA examples:
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    status, reason = is_acyl_CoA(test_smiles)
    print(status, reason)