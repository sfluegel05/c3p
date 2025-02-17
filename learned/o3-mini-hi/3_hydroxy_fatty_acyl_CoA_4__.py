"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition: An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.

This program uses substructure matching to detect two major features:
1. A 3-hydroxy acyl thioester portion, detected by either stereochemical SMARTS "[C@H](O)CC(=O)S" or "[C@@H](O)CC(=O)S".
2. A fragment characteristic of the CoA moiety (e.g. the part that contains "SCCNC(=O)CCNC(=O)") as seen in many acyl-CoA derivatives.
In addition, we require the net formal charge (computed by RDKit) to be -4.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is in the class 3-hydroxy fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string to a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First, check the net formal charge; per definition it should be -4.
    # Note: RDKit computes the formal charge based on the structure.
    if Chem.GetFormalCharge(mol) != -4:
        return False, f"Formal charge is {Chem.GetFormalCharge(mol)}, expected -4 for CoA(4-)"
    
    # Check for the presence of the 3-hydroxy acyl thioester group.
    # We define two SMARTS patterns to catch both enantiomers.
    smarts_3hydroxy_pos = "[C@H](O)CC(=O)S"
    smarts_3hydroxy_neg = "[C@@H](O)CC(=O)S"
    pattern_pos = Chem.MolFromSmarts(smarts_3hydroxy_pos)
    pattern_neg = Chem.MolFromSmarts(smarts_3hydroxy_neg)
    
    has_3hydroxy = mol.HasSubstructMatch(pattern_pos) or mol.HasSubstructMatch(pattern_neg)
    if not has_3hydroxy:
        return False, "No 3-hydroxy acyl thioester fragment ([C@H](O)CC(=O)S or [C@@H](O)CC(=O)S) found"

    # Check for a CoA fragment. In many acyl-CoA derivatives, one finds the fragment "SCCNC(=O)CCNC(=O)".
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA-specific fragment (SCCNC(=O)CCNC(=O)) found"

    # If all checks pass, we assume that the molecule is a 3-hydroxy fatty acyl-CoA(4-).
    return True, "Molecule contains a 3-hydroxy acyl thioester fragment, a CoA moiety, and has a net charge of -4"

# Example usage:
if __name__ == "__main__":
    # One of the sample SMILES strings provided (e.g. (R)-3-hydroxybutanoyl-CoA(4-))
    sample_smiles = "C[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    valid, reason = is_3_hydroxy_fatty_acyl_CoA_4__(sample_smiles)
    print(valid, reason)