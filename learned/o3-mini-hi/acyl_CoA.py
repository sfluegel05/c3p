"""
Classifies: CHEBI:17984 acyl-CoA
"""
#!/usr/bin/env python
"""
Classifies: Acyl-CoA
Definition:
  An acyl-CoA is a thioester that results from the condensation of the thiol group of coenzyme A with the carboxy group of a fatty acid.
  Its structure features an acyl thioester group (C(=O)S) directly attached to a pantetheine moiety (recognized by one of two SMARTS) that is
  in turn linked via a phosphorylated (pyrophosphate) bridge to an adenine-containing nucleotide.
Improvements:
  • Instead of simply testing for the presence of the three fragments, we now also:
      – Require that at least one thioester sulfur is either contained in or adjacent to one of the pantetheine matches.
      – Ensure that at least one pantetheine match is connected to an adenine moiety by a short (≤25 bonds) path that contains a pyrophosphate
        group (as detected by the SMARTS pattern "OP(=O)(O)OP(=O)(O)").
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    The method requires that the molecule contains:
      1. An acyl thioester fragment: "[CX3](=O)[SX2]".
      2. A pantetheine fragment which is searched for using two possible SMARTS patterns:
           a) "SCCNC(=O)CCNC(=O)"
           b) "NCCSC(=O)"
      3. An adenine moiety: "n1cnc2c(N)ncnc12"
      
    In addition, we verify that:
      A. The sulfur participating in a thioester is either directly part of a pantetheine match or at least bonded to one.
      B. At least one pantetheine match is connected to at least one adenine match by a short (≤25 bonds) path.
         The connecting path must contain a pyrophosphate group (as identified by the SMARTS pattern "OP(=O)(O)OP(=O)(O)"),
         representing the phosphate bridge found in CoA.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of (True, reason) if the molecule meets the acyl-CoA criteria,
                     otherwise (False, explanation).
    """
    # Parse the SMILES string.
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

    # 2. Locate the pantetheine fragment using two alternative SMARTS patterns.
    pantetheine_smarts_list = ["SCCNC(=O)CCNC(=O)", "NCCSC(=O)"]
    pantetheine_matches = []
    for smarts in pantetheine_smarts_list:
        pat = Chem.MolFromSmarts(smarts)
        if pat:
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

    # 4. Check that a thioester sulfur is connected to the pantetheine fragment.
    thioester_linked = False
    # In the thioester match, index 1 is the sulfur.
    for tmatch in thioester_matches:
        thio_sulfur = tmatch[1]
        for pant_match in pantetheine_matches:
            # If the thioester sulfur is directly part of the pantetheine match.
            if thio_sulfur in pant_match:
                thioester_linked = True
                break
            # Otherwise if any neighbor of the sulfur is included in a pantetheine match.
            atom = mol.GetAtomWithIdx(thio_sulfur)
            nbr_idxs = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
            if any(n in pant_match for n in nbr_idxs):
                thioester_linked = True
                break
        if thioester_linked:
            break
    if not thioester_linked:
        return False, "Thioester group is not connected to the pantetheine fragment as expected"

    # 5. Check connectivity between pantetheine and adenine fragments.
    # We require that at least one pantetheine match is connected to an adenine match by a path of length <= 25
    # which contains a pyrophosphate group. We define a pyrophosphate SMARTS pattern.
    pyrophosphate_smarts = "OP(=O)(O)OP(=O)(O)"
    pyrophosphate_pat = Chem.MolFromSmarts(pyrophosphate_smarts)
    if pyrophosphate_pat is None:
        return False, "Error creating pyrophosphate SMARTS pattern"
    
    max_path_length = 25
    found_connection = False
    # Loop over all pantetheine and adenine matches.
    for pant_match in pantetheine_matches:
        for aden_match in adenine_matches:
            for pant_atom in pant_match:
                for aden_atom in aden_match:
                    try:
                        path = Chem.GetShortestPath(mol, pant_atom, aden_atom)
                    except Exception:
                        continue
                    if path and (len(path) <= max_path_length):
                        # Extract the sub-molecule along the connecting path.
                        submol = Chem.PathToSubmol(mol, path)
                        # Check if the submol contains the pyrophosphate motif (proper phosphate bridge).
                        if submol.HasSubstructMatch(pyrophosphate_pat):
                            found_connection = True
                            break
                if found_connection:
                    break
            if found_connection:
                break
        if found_connection:
            break
    if not found_connection:
        return False, ("Pantetheine and adenine fragments are not connected via a proper pyrophosphate bridge "
                       "(a connecting path <=25 bonds lacking the 'OP(=O)(O)OP(=O)(O)' motif was found)")
    
    return True, "Contains acyl thioester attached to a pantetheine fragment and linked via a pyrophosphate bridge to an adenine moiety (acyl-CoA)"

# Example usage:
if __name__ == '__main__':
    # Test with one of the provided acyl-CoA examples.
    test_smiles = "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC([C@H](C(NCCC(NCCSC(=O)/C=C/CCCCCCCCCCCCC)=O)=O)O)(C)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"
    status, reason = is_acyl_CoA(test_smiles)
    print(status, reason)