"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition: An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups
of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.

This version revisits the prior program based on the analysis of false positives and negatives.
It:
  1) Requires that one of the two stereochemical SMARTS for the 3-hydroxy acyl thioester fragment 
     (i.e. "[C@H](O)CC(=O)S" or "[C@@H](O)CC(=O)S") is matched and the hydroxy-bearing (chiral) carbon is not in a ring.
  2) Checks that the molecule contains a CoA fragment by (a) finding a thiotemplate substructure ("SCCNC(=O)CCNC(=O)") 
     and (b) matching one of two purine (adenosine) patterns.
  3) Accepts only molecules whose computed formal charge is -4 or alternatively 0 (drawn neutral – we then add a note).
  4) Requires at least 2 oxygen atoms, each explicitly carrying –1 formal charge and directly bonded to a phosphorus atom,
     to help verify the CoA pyrophosphate unit.
     
If any check fails, a meaningful message is returned.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets all classification criteria; False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Check formal charge: accept if computed charge is -4, or (drawn as neutral) 0.
    computed_charge = Chem.GetFormalCharge(mol)
    if computed_charge not in (-4, 0):
        return False, f"Charge check failed: computed formal charge is {computed_charge} (expected -4 or 0)."

    # (2) Look for the 3-hydroxy acyl thioester fragment with the correct stereochemistry.
    # We use two SMARTS for the chiral hydroxy-bearing carbon attached to a CC(=O)S group.
    acyl_patterns = ["[C@H](O)CC(=O)S", "[C@@H](O)CC(=O)S"]
    found_valid_acyl = False
    for smarts in acyl_patterns:
        acyl_pattern = Chem.MolFromSmarts(smarts)
        if not acyl_pattern:
            continue
        matches = mol.GetSubstructMatches(acyl_pattern)
        for match in matches:
            # Ensure that the hydroxy-bearing carbon (the first atom in the match) is not in a ring.
            if not mol.GetAtomWithIdx(match[0]).IsInRing():
                found_valid_acyl = True
                break
        if found_valid_acyl:
            break
    if not found_valid_acyl:
        return False, "No valid 3-hydroxy acyl thioester fragment found (or the hydroxy-bearing carbon is in a ring)."

    # (3) Check for the CoA fragment.
    # (a) Look for the thiotemplate fragment.
    coa_thiotemplate_smarts = "SCCNC(=O)CCNC(=O)"
    coa_template = Chem.MolFromSmarts(coa_thiotemplate_smarts)
    if not mol.HasSubstructMatch(coa_template):
        return False, "No CoA-specific thiotemplate fragment (SCCNC(=O)CCNC(=O)) found."

    # (b) Look for a purine (adenosine) substructure.
    purine_patterns = ["n1cnc2ncnc12", "n1cnc2c(N)ncnc12"]
    found_purine = False
    for ps in purine_patterns:
        purine_pattern = Chem.MolFromSmarts(ps)
        if purine_pattern and mol.HasSubstructMatch(purine_pattern):
            found_purine = True
            break
    if not found_purine:
        return False, ("No adenosine/purine fragment found to support CoA identity "
                       "(tried patterns: n1cnc2ncnc12 and n1cnc2c(N)ncnc12).")

    # (4) Extra filtering: Count oxygen atoms that are directly bonded to a phosphorus and have an explicit -1 charge.
    neg_phos_oxygens = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O" and atom.GetFormalCharge() == -1:
            # Check if any neighbor is a phosphorus atom.
            if any(neighbor.GetSymbol() == "P" for neighbor in atom.GetNeighbors()):
                neg_phos_oxygens += 1
    if neg_phos_oxygens < 2:
        return False, (f"Insufficient explicit negative charges in the phosphate region "
                       f"(found {neg_phos_oxygens}; expected at least 2 to indicate the CoA 4- state).")

    # Final decision: our molecule satisfies the criteria.
    charge_note = " (NB: drawn neutral; assuming protonated form consistent with CoA(4-))" if computed_charge == 0 else ""
    return True, ("Molecule contains a valid 3-hydroxy acyl thioester moiety, a CoA fragment (thiotemplate and purine regions), "
                  "and satisfies the charge and phosphate criteria" + charge_note)

# Example usage:
if __name__ == "__main__":
    # Test with one known sample SMILES, e.g., (R)-3-hydroxybutanoyl-CoA(4-)
    sample_smiles = "C[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    valid, reason = is_3_hydroxy_fatty_acyl_CoA_4__(sample_smiles)
    print(valid, reason)