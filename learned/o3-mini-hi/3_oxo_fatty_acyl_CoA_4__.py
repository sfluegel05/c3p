"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: 3-oxo-fatty acyl-CoA(4-)

Definition:
  An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups 
  of any 3-oxo-fatty acyl-CoA.

This implementation:
  1. Confirms that the overall molecule has a net formal charge of -4.
  2. Searches for a combined motif that includes the 3-oxo fatty acyl chain (with a thioester)
     directly connected to a CoA fragment. In particular, we require the following pattern:
       "C(=O)CC(=O)SCCNC(=O)CCNC(=O)"
     which means the acyl chain motif “C(=O)CC(=O)S” is directly attached to the CoA portion,
     thus avoiding false positives where the fragments are present but not connected.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if the molecule (given by its SMILES string) belongs to the class
    3-oxo-fatty acyl-CoA(4-).

    Checks performed:
      1. The molecule must have a net formal charge of -4.
      2. The molecule must contain a combined motif that ensures its 3-oxo fatty acyl chain is 
         connected via a thioester bond to a CoA fragment. The SMARTS for this combined check is:
           "C(=O)CC(=O)SCCNC(=O)CCNC(=O)"
      3. If the combined SMARTS pattern is found, we are confident that the 3-oxo acyl part is 
         directly connected to a CoA moiety.

    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA(4-), False otherwise.
      str: Explanation for the classification decision.
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Check the overall formal charge. (According to the class definition it should be -4.)
    charge = Chem.GetFormalCharge(mol)
    if charge != -4:
        return False, f"Molecule has a formal charge of {charge}, expected -4 for CoA(4-)."

    # Define a combined SMARTS pattern that requires the 3-oxo fatty acyl motif to be directly
    # connected to a CoA fragment. 
    # Explanation of pattern:
    #  - "C(=O)"    : a carbonyl (acyl carbonyl)
    #  - "CC(=O)"   : a CH2 group followed by another carbonyl (the 3-oxo function)
    #  - "S"        : a sulfur (forming the thioester bond)
    #  - "CCNC(=O)CCNC(=O)" : a minimal heuristic for the CoA fragment starting at the sulfur.
    combined_pattern = Chem.MolFromSmarts("C(=O)CC(=O)SCCNC(=O)CCNC(=O)")
    if combined_pattern is None:
        return False, "Error constructing combined SMARTS pattern."

    # Look for the combined motif in the molecule.
    if not mol.HasSubstructMatch(combined_pattern):
        return False, "The 3-oxo fatty acyl chain motif is not properly connected to the CoA fragment as expected."

    return True, ("Molecule contains the required 3-oxo fatty acyl chain directly connected to a CoA fragment "
                  "and has a net formal charge of -4.")

# Example usage:
if __name__ == "__main__":
    # A sample SMILES from one of the provided positive examples.
    smiles_example = ("CCCC\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H]"
                      "(O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
    print("Result:", result)
    print("Reason:", reason)