"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: 3-oxo-fatty acyl-CoA(4-)
Definition:
    An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups 
    of any 3-oxo-fatty acyl-CoA.
    
This script uses heuristic substructure matching to detect both a 3-oxo fatty acyl chain motif
and a CoA moiety fragment. Additionally, it verifies that the overall molecule carries a net charge
of -4.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if the molecule (given by its SMILES) belongs to the class
    3-oxo-fatty acyl-CoA(4-).

    The algorithm checks for:
    1. A 3-oxo fatty acyl chain motif, defined here by the key substructure:
       [C(=O)]–CH2–[C(=O)]–S (SMARTS: "[CX3](=O)[CH2][CX3](=O)[SX2]").
    2. A CoA moiety fragment. In many examples the thioester carbon (S) is attached to a fragment 
       including "SCCNC(=O)CCNC(=O)". We use this fragment as a heuristic.
    3. A computed overall formal charge of −4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a 3-oxo-fatty acyl-CoA(4-), False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check for the overall formal charge (should be -4 for CoA(4-)).
    if Chem.GetFormalCharge(mol) != -4:
        return False, f"Molecule has a formal charge of {Chem.GetFormalCharge(mol)}, expected -4 for CoA(4-)."
    
    # Define a SMARTS pattern for the 3-oxo fatty acyl chain motif.
    # This pattern represents [C(=O)]–CH2–[C(=O)]–S which is key for the 3-oxo group adjacent to the thioester.
    pattern_threeoxo = Chem.MolFromSmarts("[CX3](=O)[CH2][CX3](=O)[SX2]")
    if not mol.HasSubstructMatch(pattern_threeoxo):
        return False, "Missing the key 3-oxo fatty acyl chain motif ([C(=O)]–CH2–[C(=O)]–S)."
    
    # Define a SMARTS pattern for a CoA fragment.
    # This is a heuristic pattern that picks up the part of the molecule associated with the CoA unit.
    # Many acyl-CoA SMILES include a fragment like "SCCNC(=O)CCNC(=O)".
    pattern_coa = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pattern_coa):
        return False, "Missing a CoA moiety fragment (e.g., 'SCCNC(=O)CCNC(=O)')."
    
    # If all checks pass, the molecule is classified as 3-oxo-fatty acyl-CoA(4-).
    return True, "Molecule contains the 3-oxo fatty acyl chain motif, a CoA fragment, and has a -4 charge."

# Example usage:
if __name__ == "__main__":
    # Example SMILES from one of the provided structures.
    smiles_example = "CCCC\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
    print("Result:", result)
    print("Reason:", reason)