"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: 3-oxo-fatty acyl-CoA(4-)
Definition:
    An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups 
    of any 3-oxo-fatty acyl-CoA.
    
The algorithm uses substructure searches:
  - A SMARTS for the 3-oxo fatty acyl chain motif. We use "[CX3](=O)[CX4][CX3](=O)[SX2]" 
    which should match both a methylene or a substituted linking carbon.
  - A heuristic SMARTS for a characteristic CoA fragment: "SCCNC(=O)CCNC(=O)".
  - Then, for every match of the threeoxo motif, we check if its thioester sulfur (the fourth atom)
    is included anywhere in a match for the CoA fragment. This is more robust than requiring that 
    the S atom be exactly at position 0 of the CoA match.
  - Finally, the molecule must have a net formal charge of -4.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if the molecule (given by its SMILES) belongs to the class
    3-oxo-fatty acyl-CoA(4-).

    This function checks for:
      1. A 3-oxo fatty acyl chain motif. The SMARTS "[CX3](=O)[CX4][CX3](=O)[SX2]"
         is used to allow a substituted bridging carbon.
      2. A CoA moiety fragment using a heuristic SMARTS pattern "SCCNC(=O)CCNC(=O)".
      3. That the thioester sulfur (the last atom in a threeoxo match) is contained 
         somewhere in a matching CoA fragment (ensuring connection).
      4. That the overall molecule has a formal charge of −4.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check the overall formal charge.
    charge = Chem.GetFormalCharge(mol)
    if charge != -4:
        return False, f"Molecule has a formal charge of {charge}, expected -4 for CoA(4-)."
    
    # Define and search for the 3-oxo fatty acyl chain motif.
    # This pattern covers a carbonyl (sp2 C), followed by an sp3 carbon (which may be substituted),
    # followed by another carbonyl (sp2 C) and finally the thioester sulfur.
    pattern_threeoxo = Chem.MolFromSmarts("[CX3](=O)[CX4][CX3](=O)[SX2]")
    matches_threeoxo = mol.GetSubstructMatches(pattern_threeoxo)
    if not matches_threeoxo:
        return False, "Missing the key 3-oxo fatty acyl chain motif ([C(=O)]-[CX4]-[C(=O)]-S)."
    
    # Define and search for a CoA moiety fragment.
    # The heuristic pattern here is designed to match a portion of the CoA structure.
    pattern_coa = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    matches_coa = mol.GetSubstructMatches(pattern_coa)
    if not matches_coa:
        return False, "Missing a CoA moiety fragment (e.g., 'SCCNC(=O)CCNC(=O)')."
    
    # Check that any 3-oxo motif’s thioester S atom is connected to the CoA fragment.
    # In the 3-oxo match, the fourth atom (index 3) is the S atom.
    for threeoxo_match in matches_threeoxo:
        threeoxo_sulfur = threeoxo_match[3]  # thioester sulfur from the threeoxo motif.
        # For each CoA fragment match, check if this sulfur atom is among its atoms.
        for coa_match in matches_coa:
            if threeoxo_sulfur in coa_match:
                return True, ("Molecule contains a 3-oxo fatty acyl chain motif connected "
                              "to the CoA fragment and has a net -4 charge.")
    
    # If none of the threeoxo sulfur atoms overlap with any CoA fragment, then the motifs are not connected.
    return False, "The 3-oxo fatty acyl chain motif is not connected to the CoA fragment as expected."

# Example usage:
if __name__ == "__main__":
    # One example SMILES from the provided positive list.
    smiles_example = ("CCCC\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H]"
                      "(O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
    print("Result:", result)
    print("Reason:", reason)