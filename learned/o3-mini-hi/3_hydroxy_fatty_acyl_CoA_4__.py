"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition: An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups
of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.

This program uses substructure matching to detect three key features:
1. A 3-hydroxy acyl thioester fragment (in either stereochemical form).
2. A CoA moiety: we require (a) a CoA “thiotemplate” fragment and (b) a purine (adenosine) fragment.
   Here, we now check for either an unsubstituted purine ring pattern or one with a substitution.
3. Charge: either the overall molecular formal charge equals –4 OR at least 4 negatively charged oxygens
   bound to phosphorus are detected.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule belongs to the 3-hydroxy fatty acyl-CoA(4-) class.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----------------------------
    # Charge criteria:
    # Either the overall computed formal charge is –4 OR
    # at least 4 oxygens with -1 charge directly bound to a phosphorus are found.
    computed_charge = Chem.GetFormalCharge(mol)
    charge_ok = False
    if computed_charge == -4:
        charge_ok = True
    else:
        count_neg_o_on_p = 0
        for atom in mol.GetAtoms():
            # Look for oxygen with -1 formal charge
            if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1:
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 15:  # phosphorus
                        count_neg_o_on_p += 1
                        break
        if count_neg_o_on_p >= 4:
            charge_ok = True
    if not charge_ok:
        return False, f"Charge check failed: computed formal charge is {computed_charge} and found only {(count_neg_o_on_p if computed_charge != -4 else 'N/A')} negatively charged oxygen(s) on phosphorus; expected -4 overall or at least 4 such oxygens."

    # ----------------------------
    # Check for the 3-hydroxy acyl thioester fragment.
    # Two SMARTS patterns for the two stereochemical variants.
    smarts_3hydroxy_pos = "[C@H](O)CC(=O)S"
    smarts_3hydroxy_neg = "[C@@H](O)CC(=O)S"
    pattern_pos = Chem.MolFromSmarts(smarts_3hydroxy_pos)
    pattern_neg = Chem.MolFromSmarts(smarts_3hydroxy_neg)
    if not (mol.HasSubstructMatch(pattern_pos) or mol.HasSubstructMatch(pattern_neg)):
        return False, "No 3-hydroxy acyl thioester fragment ([C@H](O)CC(=O)S or [C@@H](O)CC(=O)S) found."
    
    # ----------------------------
    # Check for the CoA fragment.
    # (a) Look for the CoA thiotemplate fragment.
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA-specific thiotemplate fragment (SCCNC(=O)CCNC(=O)) found."
    
    # (b) Look for an adenosine/purine fragment.
    # Use two SMARTS patterns:
    # Pattern 1: unsubstituted purine ring as originally specified.
    purine_smarts1 = "n1cnc2ncnc12"
    # Pattern 2: purine ring with a substitution on the carbon (e.g., n1cnc2c(N)ncnc12)
    purine_smarts2 = "n1cnc2c(N)ncnc12"
    purine_pattern1 = Chem.MolFromSmarts(purine_smarts1)
    purine_pattern2 = Chem.MolFromSmarts(purine_smarts2)
    if not (mol.HasSubstructMatch(purine_pattern1) or mol.HasSubstructMatch(purine_pattern2)):
        return False, "No adenosine/purine fragment found to support CoA identity (tried patterns: n1cnc2ncnc12 and n1cnc2c(N)ncnc12)."
    
    # All tests passed.
    return True, "Molecule contains a 3-hydroxy acyl thioester moiety, a CoA moiety (both thiotemplate and purine fragments), and meets the charge criteria."

# Example usage:
if __name__ == "__main__":
    # Test with one sample SMILES, e.g., for (R)-3-hydroxybutanoyl-CoA(4-)
    sample_smiles = "C[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    valid, reason = is_3_hydroxy_fatty_acyl_CoA_4__(sample_smiles)
    print(valid, reason)