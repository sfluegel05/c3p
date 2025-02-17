"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition: An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups
of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.

This improved program uses substructure matching for the 3-hydroxy acyl thioester fragment (ensuring that 
the hydroxy-bearing chiral carbon is not in a ring), verifies the presence of key CoA fragments (both 
a thiotemplate and an adenosine/purine part), checks that the formal charge is either –4 or 0 (to allow for 
neutral-drawn SMILES), and additionally requires that there are at least two oxygen atoms explicitly carrying 
a -1 charge in the phosphate/pyrophosphate region (to avoid false positives arising from neutral phosphate groups).
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Our revised procedure relies on:
      1. Detecting the 3-hydroxy acyl thioester fragment. We look for one (or both) of the patterns:
         "[C@H](O)CC(=O)S" or "[C@@H](O)CC(=O)S". In addition, we require that the hydroxy-bearing 
         (chiral) carbon is not part of any ring.
      2. Confirming the presence of a CoA fragment by:
         (a) Finding a “thiotemplate” substructure: "SCCNC(=O)CCNC(=O)".
         (b) Detecting a purine (adenosine) fragment via one of two patterns: 
             "n1cnc2ncnc12" or "n1cnc2c(N)ncnc12".
      3. Applying a charge check: while many are drawn as neutral (formal charge 0), the actual species 
         should have a total charge of -4. We accept if the computed formal charge is either -4 or 0.
      4. Added filtering: We count the oxygen atoms that are directly bound to phosphorus and have an explicit 
         formal charge of –1. In our true CoA species (with the expected pyrophosphate unit), we require at least 2.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the classification.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ----------------------------
    # (1) Charge criteria:
    # Although many structures are drawn in their neutral form, we accept only molecules whose computed formal charge 
    # is either -4 or 0.
    computed_charge = Chem.GetFormalCharge(mol)
    if computed_charge not in (-4, 0):
        return False, f"Charge check failed: computed formal charge is {computed_charge} (expected -4 or 0)"
    
    # ----------------------------
    # (2) Check for the 3-hydroxy acyl thioester fragment.
    # Look for either stereochemical version and then ensure that the hydroxy-bearing carbon is not in a ring.
    acyl_smarts = ["[C@H](O)CC(=O)S", "[C@@H](O)CC(=O)S"]
    found_valid_acyl = False
    for smarts in acyl_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            chiral_atom = mol.GetAtomWithIdx(match[0])
            if not chiral_atom.IsInRing():
                found_valid_acyl = True
                break
        if found_valid_acyl:
            break
    if not found_valid_acyl:
        return False, "No valid 3-hydroxy acyl thioester fragment found (or the hydroxy-bearing carbon is in a ring)."
    
    # ----------------------------
    # (3) Check for the CoA fragment.
    # (a) Check for the thiotemplate fragment.
    coa_thiotemplate_smarts = "SCCNC(=O)CCNC(=O)"
    coa_template = Chem.MolFromSmarts(coa_thiotemplate_smarts)
    if not mol.HasSubstructMatch(coa_template):
        return False, "No CoA-specific thiotemplate fragment (SCCNC(=O)CCNC(=O)) found."
    
    # (b) Check for a purine part (adenosine) within the molecule.
    purine_smarts = ["n1cnc2ncnc12", "n1cnc2c(N)ncnc12"]
    found_purine = False
    for ps in purine_smarts:
        purine_pattern = Chem.MolFromSmarts(ps)
        if purine_pattern and mol.HasSubstructMatch(purine_pattern):
            found_purine = True
            break
    if not found_purine:
        return False, ("No adenosine/purine fragment found to support CoA identity "
                       "(tried patterns: n1cnc2ncnc12 and n1cnc2c(N)ncnc12).")
    
    # ----------------------------
    # (4) Extra filtering on the phosphate groups.
    # We inspect all atoms and count oxygen atoms that are directly bonded to a phosphorus and carry a -1 formal charge.
    neg_phosphate_oxygens = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and atom.GetFormalCharge() == -1:
            # Check if any neighbor is phosphorus:
            if any(neigh.GetSymbol() == "P" for neigh in atom.GetNeighbors()):
                neg_phosphate_oxygens += 1
    if neg_phosphate_oxygens < 2:
        return False, ("Insufficient explicit negative charges in the phosphate region "
                       f"(found {neg_phosphate_oxygens}; expected at least 2 to indicate the CoA 4- state).")
    
    # ----------------------------
    # Final decision.
    charge_note = " (NB: drawn neutral; assuming protonated form consistent with CoA(4-))" if computed_charge == 0 else ""
    return True, ("Molecule contains a valid 3-hydroxy acyl thioester moiety, a CoA fragment (thiotemplate and purine regions), "
                  "and satisfies the charge and phosphate criteria" + charge_note)

# Example usage:
if __name__ == "__main__":
    # Test with one known sample SMILES, e.g., (R)-3-hydroxybutanoyl-CoA(4-)
    sample_smiles = "C[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    valid, reason = is_3_hydroxy_fatty_acyl_CoA_4__(sample_smiles)
    print(valid, reason)