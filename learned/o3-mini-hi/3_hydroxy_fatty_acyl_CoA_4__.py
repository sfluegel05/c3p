"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition: An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups
of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.

This updated program uses substructure matching to detect three key features:
1. A 3-hydroxy acyl thioester fragment. (Both stereochemical forms are checked via SMARTS.)
2. A CoA moiety: we now require both the “thiotemplate” fragment and the detection of an adenosine/purine fragment.
3. Charge: Either the overall computed formal charge equals –4 OR we can detect at least 4 negatively charged oxygens attached to phosphorus.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is in the class 3-hydroxy fatty acyl-CoA(4-), False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----------------------------
    # Check the “charge” condition.
    # Either the overall computed formal charge should be –4
    # OR we count at least 4 oxygen atoms with explicit charge –1 that are bound to a phosphorus.
    computed_charge = Chem.GetFormalCharge(mol)
    charge_ok = False
    if computed_charge == -4:
        charge_ok = True
    else:
        # Count oxygen atoms with formal charge -1 that are directly connected to a phosphorus atom.
        count_neg_o_on_p = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 15:  # phosphorus
                        count_neg_o_on_p += 1
                        break
        if count_neg_o_on_p >= 4:
            charge_ok = True
    if not charge_ok:
        return False, f"Charge check failed: computed formal charge is {computed_charge} and found only {count_neg_o_on_p if computed_charge != -4 else 'N/A'} negatively charged oxygen(s) on phosphorus; expected -4 overall in a proper CoA(4-)"
    
    # ----------------------------
    # Check for the 3-hydroxy acyl thioester fragment.
    # Use two SMARTS strings for the two stereochemical variants.
    smarts_3hydroxy_pos = "[C@H](O)CC(=O)S"
    smarts_3hydroxy_neg = "[C@@H](O)CC(=O)S"
    pattern_pos = Chem.MolFromSmarts(smarts_3hydroxy_pos)
    pattern_neg = Chem.MolFromSmarts(smarts_3hydroxy_neg)
    if not (mol.HasSubstructMatch(pattern_pos) or mol.HasSubstructMatch(pattern_neg)):
        return False, "No 3-hydroxy acyl thioester fragment ([C@H](O)CC(=O)S or [C@@H](O)CC(=O)S) found"
    
    # ----------------------------
    # Check for the CoA fragment.
    # We require two aspects:
    # a) a fragment containing the thiotemplate part: "SCCNC(=O)CCNC(=O)"
    # b) an adenosine (purine) fragment: common acyl-CoA SMILES include a fragment like "n1cnc2ncnc12"
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    purine_smarts = "n1cnc2ncnc12"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    purine_pattern = Chem.MolFromSmarts(purine_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA-specific thiotemplate fragment (SCCNC(=O)CCNC(=O)) found"
    if not mol.HasSubstructMatch(purine_pattern):
        return False, "No adenosine/purine fragment (n1cnc2ncnc12) found to support CoA identity"
    
    # If all tests pass, the molecule is most likely a 3-hydroxy fatty acyl-CoA(4-)
    return True, "Molecule contains a 3-hydroxy acyl thioester moiety, a CoA moiety (both thiotemplate and purine fragment), and the charge criteria (–4 overall or sufficient negative oxygen atoms) are met."

# Example usage:
if __name__ == "__main__":
    # Test with one sample SMILES:
    sample_smiles = "C[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    valid, reason = is_3_hydroxy_fatty_acyl_CoA_4__(sample_smiles)
    print(valid, reason)