"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: 3-oxo-fatty acyl-CoA(4-)
Definition:
    An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups 
    of any 3-oxo-fatty acyl-CoA.
    
Improvements over the previous attempt:
  - The SMARTS pattern for the 3-oxo fatty acyl chain motif now allows for substitution 
    on the central carbon (using [CX4] instead of [CH2]).
  - We retrieve all substructure matches for both the 3-oxo motif and the CoA fragment.
  - We then check that the thioester sulfur (last atom in the 3-oxo motif match) overlaps 
    with the sulfur of the CoA fragment (first atom in its match). This ensures that the 3-oxo 
    motif is directly connected to the CoA moiety.
  - The overall molecule must have a formal charge of -4.
    
A correct classification now means the molecule contains a 3-oxo motif that is connected 
to the CoA fragment and has the correct net charge.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if the molecule (given by its SMILES) belongs to the class
    3-oxo-fatty acyl-CoA(4-).

    The algorithm checks for:
      1. A 3-oxo fatty acyl chain motif. Previously, a rigid pattern "[CX3](=O)[CH2][CX3](=O)[SX2]" 
         was used. Here, we allow substitution on the central carbon by using "[CX4]" so that 
         both CH2 and substituted CH groups match.
      2. A CoA moiety fragment using a heuristic pattern "SCCNC(=O)CCNC(=O)".
      3. Overlap of the sulfur atom in the two matches. In the 3-oxo motif the fourth atom 
         corresponds to the thioester sulfur ([SX2]), which should be the same as the first atom 
         found in the CoA fragment match.
      4. A computed overall formal charge of −4.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Verify that the overall formal charge is -4.
    charge = Chem.GetFormalCharge(mol)
    if charge != -4:
        return False, f"Molecule has a formal charge of {charge}, expected -4 for CoA(4-)."

    # Define a SMARTS pattern for the 3-oxo fatty acyl chain motif.
    # Previously: "[CX3](=O)[CH2][CX3](=O)[SX2]"
    # Updated to allow substitution on the linking carbon: "[CX3](=O)[CX4][CX3](=O)[SX2]"
    pattern_threeoxo = Chem.MolFromSmarts("[CX3](=O)[CX4][CX3](=O)[SX2]")
    matches_threeoxo = mol.GetSubstructMatches(pattern_threeoxo)
    if not matches_threeoxo:
        return False, "Missing the key 3-oxo fatty acyl chain motif ([C(=O)]–[CX4]–[C(=O)]–S)."
    
    # Define a SMARTS pattern for a CoA fragment.
    # This heuristic pattern captures a common fragment of CoA: "SCCNC(=O)CCNC(=O)"
    pattern_coa = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    matches_coa = mol.GetSubstructMatches(pattern_coa)
    if not matches_coa:
        return False, "Missing a CoA moiety fragment (e.g., 'SCCNC(=O)CCNC(=O)')."
    
    # Check for an overlap between the two motifs.
    # For pattern_threeoxo the order is: [C(=O)] [CX4] [C(=O)] [S]
    # Hence, the thioester S atom is at index 3 of the match tuple.
    # For pattern_coa the first atom (index 0) is the S atom.
    for to_match in matches_threeoxo:
        threeoxo_S = to_match[3]
        for coa_match in matches_coa:
            coa_S = coa_match[0]
            if threeoxo_S == coa_S:
                return True, ("Molecule contains the 3-oxo fatty acyl chain motif connected to "
                              "the CoA fragment and has a -4 charge.")
    
    # If no overlapping S atom is found, the motifs are not properly connected.
    return False, "The 3-oxo fatty acyl chain motif is not connected to the CoA fragment as expected."

# Example usage:
if __name__ == "__main__":
    # One example SMILES from the positive list.
    smiles_example = ("CCCC\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H]"
                      "(O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
    print("Result:", result)
    print("Reason:", reason)