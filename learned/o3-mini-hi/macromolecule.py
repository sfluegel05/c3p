"""
Classifies: CHEBI:33839 macromolecule
"""
#!/usr/bin/env python3
"""
Classifies: Macromolecule

Definition:
  "A macromolecule is a molecule of high relative molecular mass, 
   the structure of which essentially comprises the multiple repetition 
   of units derived, actually or conceptually, from molecules of low relative molecular mass."

Improvement:
  • We observed that simply checking for high weight and minimal “repetition evidence”
    (i.e. one or two amide bonds or sugar motifs) caused many false positives.
  • We now require that the sum of detected repeating units (amide bonds + sugar‐like rings)
    meets stronger minimum thresholds that vary with molecular weight range and overall atom count.
  • The thresholds are as follows:
         - If weight ≥ 1000 Da, require at least 2 repeating units AND at least 40 atoms.
         - For weight between 800 and 1000 Da, require at least 3 repeating units AND at least 50 atoms.
         - For weight between 500 and 800 Da, require at least 4 repeating units AND at least 60 atoms.
         - Otherwise, classify as not macromolecular.
      
Usage: Call is_macromolecule(smiles) to get a boolean plus a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    
    Revisions based on feedback:
      - Incorporate a scaled threshold on repetition evidence (sum of amide bonds and sugar units)
        plus a minimum overall atom count.
        
      Heuristic:
         • If molecular weight ≥ 1000 Da:
              require (amide bonds + sugar units) >= 2  AND  total atoms >= 40.
         • If weight is in [800, 1000) Da:
              require (amide bonds + sugar units) >= 3  AND  total atoms >= 50.
         • If weight is in [500, 800) Da:
              require (amide bonds + sugar units) >= 4  AND  total atoms >= 60.
         • Otherwise, the molecule is not classified as a macromolecule.
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        (bool, str): Classification result and explanation.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight and total atom count.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    num_atoms = mol.GetNumAtoms()
    
    # Count amide bonds.
    # A simple SMARTS "C(=O)N" can detect many peptide-type bonds.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    
    # Detect sugar-like units using a pyranose-like pattern.
    # (A more refined approach would allow for furanoses and different stereochemistries.)
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1OC(O)C(O)C(O)C1O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    num_sugar = len(sugar_matches)
    
    # Combine repetition evidence.
    repetition_evidence = num_amide + num_sugar
    
    # Apply the heuristic thresholds.
    if mol_wt >= 1000:
        if repetition_evidence >= 2 and num_atoms >= 40:
            reason = (f"Molecular weight is {mol_wt:.1f} Da which is high for a macromolecule. "
                      f"Found {num_amide} amide bonds and {num_sugar} sugar units (total evidence: {repetition_evidence}), "
                      f"with {num_atoms} atoms overall.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da but repetition evidence is low "
                      f"(only {num_amide} amide bonds + {num_sugar} sugar units = {repetition_evidence}) "
                      f"and atom count is {num_atoms}; criteria for repeated subunits are not met.")
            return False, reason
    elif 800 <= mol_wt < 1000:
        if repetition_evidence >= 3 and num_atoms >= 50:
            reason = (f"Molecular weight is {mol_wt:.1f} Da, and the molecule has {num_atoms} atoms "
                      f"with repetition evidence of {repetition_evidence} (from {num_amide} amide bonds and {num_sugar} sugar units), "
                      "suggesting repeated units.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da but repetition evidence is only {repetition_evidence} "
                      f"and atom count is {num_atoms}; criteria for repeated subunits are not met.")
            return False, reason
    elif 500 <= mol_wt < 800:
        if repetition_evidence >= 4 and num_atoms >= 60:
            reason = (f"Molecular weight is {mol_wt:.1f} Da, and the molecule has {num_atoms} atoms "
                      f"with repetition evidence of {repetition_evidence} (from {num_amide} amide bonds and {num_sugar} sugar units), "
                      "suggesting multiple repeat units.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da with {num_atoms} atoms and only {repetition_evidence} "
                      f"(amide bonds + sugar units) found; criteria for a macromolecule are not met.")
            return False, reason
    else:
        reason = (f"Molecular weight is {mol_wt:.1f} Da with {num_atoms} atoms and "
                  f"{repetition_evidence} repeating units (from {num_amide} amide bonds and {num_sugar} sugar units); "
                  "criteria for a macromolecule are not met.")
        return False, reason

# Example usage:
if __name__ == '__main__':
    # Test with a known peptide macromolecule SMILES.
    test_smiles = "CC[C@H](C)[C@H]1C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)NCCCC[C@@H](C(=O)N[C@@H](C(=O)N1)CCCN)NC(=O)[C@H]([C@@H](C)CC)NC(=O)[C@@H](CCC(=O)O)NC(=O)[C@H](CC(C)C)NC(=O)C2CSC(=N2)C(C(C)CC)N)CC(=O)N)CC(=O)O)CC3=CN=CN3)CC4=CC=CC=C4"
    flag, explanation = is_macromolecule(test_smiles)
    print("Is macromolecule?", flag)
    print("Reason:", explanation)