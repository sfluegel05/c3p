"""
Classifies: CHEBI:15489 3-oxo-fatty acyl-CoA
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-fatty acyl-CoA
An oxo fatty acyl-CoA results from the condensation of the thiol group of Coenzyme A
with the carboxy group of any 3-oxo fatty acid.
This program searches for a 3-oxo fatty acyl thioester fragment. Two alternate SMARTS are used:
   1. R–C(=O)–CH2–C(=O)S
   2. S–C(=O)–CH2–C(=O)–R
in order to capture both orientations of the thioester.
Additionally, a Coenzyme A adenine fragment (n1cnc2ncnc(c12)) must be present.
"""

from rdkit import Chem

def is_3_oxo_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA based on its SMILES string.
    
    Criteria:
      1. The SMILES must parse and be neutral (no deprotonated oxygens "[O-]").
      2. The molecule must contain a Coenzyme A adenine fragment.
      3. The molecule must contain a 3-oxo fatty acyl thioester fragment.
         We look for one of two possible fragments:
             a) R–C(=O)–CH2–C(=O)S
             b) S–C(=O)–CH2–C(=O)–R
         This fragment corresponds to the 3-oxo group (at the 3‐position of a fatty acid)
         and the thioester bond linking to CoA.
    
    Args:
      smiles (str): SMILES string representing the molecule.
      
    Returns:
      bool: True if the molecule meets the criteria, False otherwise.
      str: Explanation for the classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # We require neutral compounds (no deprotonated oxygens).
    if "[O-]" in smiles:
        return False, "Contains deprotonated oxygens; expected a neutral CoA derivative"
    
    # Check for a Coenzyme A adenine fragment.
    coa_pattern = Chem.MolFromSmarts("n1cnc2ncnc(c12)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Does not contain a Coenzyme A adenine fragment"
    
    # Define two SMARTS patterns for the 3-oxo fatty acyl thioester fragment.
    # Pattern 1: R–C(=O)–CH2–C(=O)S  (e.g. acetoacetyl-CoA: CH3-C(=O)-CH2-C(=O)S)
    frag_smarts1 = "[#6]-C(=O)[CH2]C(=O)S"
    # Pattern 2: S–C(=O)–CH2–C(=O)–R  (if the orientation is reversed)
    frag_smarts2 = "S-C(=O)[CH2]C(=O)-[#6]"
    
    frag_patterns = [frag_smarts1, frag_smarts2]
    valid_fragment_found = False
    reason_details = ""
    
    for smarts in frag_patterns:
        frag_pattern = Chem.MolFromSmarts(smarts)
        if frag_pattern is None:
            continue  # should not happen, but safety first
        matches = mol.GetSubstructMatches(frag_pattern)
        if matches:
            valid_fragment_found = True
            break  # At least one valid fragment is found.
    
    if not valid_fragment_found:
        reason_details = "Does not contain an uncompromised 3-oxo fatty acyl thioester substructure"
        return False, reason_details
    
    return True, "Contains a valid 3-oxo fatty acyl thioester fragment linked to a Coenzyme A moiety"

# Example usage:
if __name__ == "__main__":
    # Test using one of the provided sample SMILES: 3-oxoadipyl-CoA
    example_smiles = ("CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)"
                        "[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)CC(=O)CCC(O)=O")
    result, reason = is_3_oxo_fatty_acyl_CoA(example_smiles)
    print(result, reason)