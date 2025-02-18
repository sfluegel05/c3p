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
  2. Searches for a 3-oxo fatty acyl chain motif with a thioester linkage. For this, 
     we use a SMARTS pattern "C(=O)CC(=O)S" that covers an acyl carbonyl attached to a 
     methylene (CH2), followed by another carbonyl and then the sulfur.
  3. Searches for a CoA fragment using a heuristic SMARTS pattern "SCCNC(=O)CCNC(=O)".
  4. Verifies that the thioester sulfur (atom index 3 in the acyl chain motif) is 
     contained within one of the CoA fragment matches.

If all these conditions pass the function returns True plus an explanation; otherwise,
it returns False with a reason.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if the molecule (given by its SMILES) belongs to the class
    3-oxo-fatty acyl-CoA(4-).

    Checks performed:
      1. The molecule must have a net formal charge of -4.
      2. It must contain a 3-oxo fatty acyl chain motif connected via a thioester bond,
         recognized by the SMARTS "C(=O)CC(=O)S".
      3. It must contain a CoA moiety fragment, here approximated by the SMARTS
         "SCCNC(=O)CCNC(=O)".
      4. The sulfur atom from the acyl chain motif (the thioester S) must be part
         of the CoA fragment, ensuring the two parts are connected.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA(4-), False otherwise.
      str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Verify the overall formal charge.
    charge = Chem.GetFormalCharge(mol)
    if charge != -4:
        return False, f"Molecule has a formal charge of {charge}, expected -4 for CoA(4-)."

    # Define the 3-oxo fatty acyl chain motif with thioester linkage.
    # This pattern represents: carbonyl (C(=O)) - CH2 - carbonyl (C(=O)) - sulfur (S)
    pattern_chain = Chem.MolFromSmarts("C(=O)CC(=O)S")
    matches_chain = mol.GetSubstructMatches(pattern_chain)
    if not matches_chain:
        return False, "Missing the 3-oxo fatty acyl chain motif with thioester linkage (C(=O)CC(=O)S)."

    # Define a heuristic fragment for the CoA moiety.
    pattern_coa = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    matches_coa = mol.GetSubstructMatches(pattern_coa)
    if not matches_coa:
        return False, "Missing a CoA moiety fragment (e.g., 'SCCNC(=O)CCNC(=O)')."

    # For each occurrence of the 3-oxo motif, check that its thioester sulfur is connected to the CoA fragment.
    # In our pattern "C(=O)CC(=O)S", the sulfur is the fourth atom (index 3).
    for chain_match in matches_chain:
        chain_sulfur = chain_match[3]
        # Verify if the sulfur atom is part of any CoA fragment match.
        for coa_match in matches_coa:
            if chain_sulfur in coa_match:
                return True, ("Molecule contains a 3-oxo fatty acyl chain motif connected to the CoA fragment "
                              "and has a net formal charge of -4.")
    
    return False, "The 3-oxo fatty acyl chain motif is not connected to the CoA fragment as expected."

# Example usage:
if __name__ == "__main__":
    # One example SMILES from the provided list of positive cases:
    smiles_example = ("CCCC\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H]"
                      "(O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
    print("Result:", result)
    print("Reason:", reason)