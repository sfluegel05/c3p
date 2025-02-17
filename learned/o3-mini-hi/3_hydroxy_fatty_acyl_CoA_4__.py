"""
Classifies: CHEBI:65102 3-hydroxy fatty acyl-CoA(4-)
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA(4-)
Definition: An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups
of any 3-hydroxy fatty acyl-CoA; major species at pH 7.3.

This revised program uses substructure matching and a relaxed charge test (to allow for neutral SMILES)
plus extra filtering on the 3-hydroxy fragment to avoid cases where the hydroxyl is embedded in a ring.
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA(4-) based on its SMILES string.
    
    Our classification procedure relies on:
      1. Detecting the 3-hydroxy acyl thioester fragment. We look for one (or both) of the patterns:
         "[C@H](O)CC(=O)S" or "[C@@H](O)CC(=O)S". In addition, we filter out any match where the chiral,
         hydroxy-bearing carbon is in a ring (to avoid ring-bound hydroxyls seen in pyran fragments).
      2. Confirming the presence of a CoA fragment by:
         (a) Finding a “thiotemplate” substructure: "SCCNC(=O)CCNC(=O)".
         (b) Detecting a purine (adenosine) fragment via one of two patterns: 
             "n1cnc2ncnc12" or "n1cnc2c(N)ncnc12".
      3. Applying a charge check. Many of these molecules commonly appear drawn as neutral (formal charge 0)
         even though they are the 4- species; therefore we accept molecules when the computed formal charge is
         either exactly –4 or 0. For any other value the molecule is rejected.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the classification.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # -------------------------------------
    # Charge criteria:
    # In many cases the SMILES are drawn in their neutral form.
    # We accept if the computed formal charge is either -4 or 0.
    computed_charge = Chem.GetFormalCharge(mol)
    if computed_charge not in (-4, 0):
        return False, f"Charge check failed: computed formal charge is {computed_charge} (expected -4 or 0 for protonated form)."
    
    # -------------------------------------
    # Check for the 3-hydroxy acyl thioester fragment.
    # Use two SMARTS patterns that differ by stereochemistry.
    smarts_patterns = ["[C@H](O)CC(=O)S", "[C@@H](O)CC(=O)S"]
    found_valid_acyl = False
    for smarts in smarts_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            # The first atom in our pattern is the chiral (hydroxy-bearing) carbon.
            # We require that this atom is not part of any ring.
            chiral_atom = mol.GetAtomWithIdx(match[0])
            if not chiral_atom.IsInRing():
                found_valid_acyl = True
                break
        if found_valid_acyl:
            break
    if not found_valid_acyl:
        return False, "No valid 3-hydroxy acyl thioester fragment found (or the hydroxy-bearing carbon is in a ring)."
    
    # -------------------------------------
    # Check for the CoA fragment.
    # (a) CoA thiotemplate fragment (allows some variability).
    coa_thiotemplate_smarts = "SCCNC(=O)CCNC(=O)"
    coa_template = Chem.MolFromSmarts(coa_thiotemplate_smarts)
    if not mol.HasSubstructMatch(coa_template):
        return False, "No CoA-specific thiotemplate fragment (SCCNC(=O)CCNC(=O)) found."
    
    # (b) Purine (adenosine) fragment. Two variants are allowed.
    purine_smarts = ["n1cnc2ncnc12", "n1cnc2c(N)ncnc12"]
    found_purine = False
    for ps in purine_smarts:
        purine_pattern = Chem.MolFromSmarts(ps)
        if purine_pattern and mol.HasSubstructMatch(purine_pattern):
            found_purine = True
            break
    if not found_purine:
        return False, "No adenosine/purine fragment found to support CoA identity (tried patterns: n1cnc2ncnc12 and n1cnc2c(N)ncnc12)."
    
    # -------------------------------------
    # Final conclusion.
    # If the formal charge was 0 (i.e. drawn in a protonated form) there is a note.
    charge_note = " (NB: drawn neutral; assuming protonated form consistent with CoA(4-))" if computed_charge == 0 else ""
    return True, f"Molecule contains a valid 3-hydroxy acyl thioester moiety, a CoA fragment (thiotemplate and purine regions), and satisfies the charge criteria{charge_note}."

# Example usage:
if __name__ == "__main__":
    # Test with one known sample SMILES, e.g., for (R)-3-hydroxybutanoyl-CoA(4-)
    sample_smiles = "C[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    valid, reason = is_3_hydroxy_fatty_acyl_CoA_4__(sample_smiles)
    print(valid, reason)