"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition:
  An acyl-CoA is a thioester that results from condensation of the thiol group of Coenzyme A
  with the carboxy group of a fatty acid. Its structure features an acyl thioester group (C(=O)S)
  directly attached to a pantetheine moiety (here recognized by a fragment “SCCNC(=O)CCNC(=O)”)
  that in turn is connected (via a phosphate bridge) to an adenine-containing nucleotide.
Improvement:
  Instead of only testing for the mere presence of two fragments we:
    1. Look for an acyl thioester group "[CX3](=O)[SX2]".
    2. Look for a pantetheine fragment "SCCNC(=O)CCNC(=O)".
    3. Look for an adenine moiety "n1cnc2c(N)ncnc12".
    4. Check that the thioester sulfur is part of (or immediately adjacent to) the pantetheine match.
    5. Check that the pantetheine fragment and the adenine moiety are connected via a short bond‐path (<=25 bonds)
       that contains at least one phosphorus (P) atom (i.e. a bridging phosphate).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    The molecule must contain:
      1. An acyl thioester fragment: "[CX3](=O)[SX2]".
      2. A pantetheine fragment: "SCCNC(=O)CCNC(=O)".
      3. An adenine moiety: "n1cnc2c(N)ncnc12".
    In addition, we verify that:
      A. The thioester sulfur (from the acyl fragment) is part of -- or directly adjacent to --
         a pantetheine match.
      B. At least one pantetheine match is connected to the adenine moiety by a short topological
         distance (<=25 bonds) and that the connecting atoms include at least one phosphorus (P) atom.
    
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
    
    # 1. Match the acyl thioester group. (It must appear as a carbonyl attached to sulfur.)
    thioester_smarts = "[CX3](=O)[SX2]"
    thioester_pat = Chem.MolFromSmarts(thioester_smarts)
    if thioester_pat is None:
        return False, "Error creating thioester SMARTS"
    thioester_matches = mol.GetSubstructMatches(thioester_pat)
    if not thioester_matches:
        return False, "Thioester group '[CX3](=O)[SX2]' not found"
    
    # 2. Match the pantetheine fragment.
    pantetheine_smarts = "SCCNC(=O)CCNC(=O)"
    pantetheine_pat = Chem.MolFromSmarts(pantetheine_smarts)
    if pantetheine_pat is None:
        return False, "Error creating pantetheine SMARTS pattern"
    pantetheine_matches = mol.GetSubstructMatches(pantetheine_pat)
    if not pantetheine_matches:
        return False, "Pantetheine fragment 'SCCNC(=O)CCNC(=O)' not detected"
    
    # 3. Match the adenine moiety.
    adenine_smarts = "n1cnc2c(N)ncnc12"
    adenine_pat = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pat is None:
        return False, "Error creating adenine SMARTS pattern"
    adenine_matches = mol.GetSubstructMatches(adenine_pat)
    if not adenine_matches:
        return False, "Adenine moiety of CoA not detected"
    
    # 4. Check that a thioester group is connected to the pantetheine.
    # We require that the sulfur atom in an acyl thioester match is present within or immediately adjacent (neighbors)
    # to one of the atoms in a pantetheine match.
    thioester_linked = False
    for match in thioester_matches:
        # In our SMARTS "[CX3](=O)[SX2]", index 0 is the carbonyl C, index 1 is the sulfur.
        thioester_sulfur = match[1]
        for pant_match in pantetheine_matches:
            # If the thioester sulfur atom is in the pantetheine match, consider it connected
            if thioester_sulfur in pant_match:
                thioester_linked = True
                break
            # Alternatively, if it has a direct bond to an atom in the pantetheine match.
            nbrs = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(thioester_sulfur).GetNeighbors()]
            if any(n in pant_match for n in nbrs):
                thioester_linked = True
                break
        if thioester_linked:
            break
    if not thioester_linked:
        return False, "Thioester group is not connected to the pantetheine fragment as expected"
    
    # 5. Check connectivity between the pantetheine fragment and the adenine moiety.
    # For each combination, compute the shortest topological path.
    # We also require that the connecting path is not too long (<= 25 bonds)
    # and that the path contains at least one phosphorus atom (representing a bridging phosphate group).
    max_path_length = 25
    found_connection = False
    for pant_match in pantetheine_matches:
        for adenine_match in adenine_matches:
            # Look at each atom pair (one from pantetheine, one from adenine)
            for pant_atom in pant_match:
                for aden_atom in adenine_match:
                    try:
                        path = Chem.GetShortestPath(mol, pant_atom, aden_atom)
                    except Exception:
                        continue
                    if path and len(path) <= max_path_length:
                        # Check if at least one atom along the path is a phosphorus (atomic number 15).
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
        return False, "Pantetheine and adenine fragments are not connected via a proper phosphate linker"
    
    # Passed all tests.
    return True, "Contains acyl thioester attached to a pantetheine fragment and linked to an adenine moiety (acyl-CoA)"

# Example usage:
# status, reason = is_acyl_CoA("SMILES_STRING_HERE")
# print(status, reason)

# For testing one of the true positives:
if __name__ == '__main__':
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    status, reason = is_acyl_CoA(test_smiles)
    print(status, reason)