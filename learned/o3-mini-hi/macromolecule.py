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

This version uses several features:
  • Molecular weight (using RDKit’s exact molecular weight).
  • Total atom count as a proxy for overall size.
  • Number of amide bonds (using a SMARTS search for the C(=O)N motif).
  • Number of repeated sugar-like units (detected via a common hexose pattern).
Based on these, different conditions are applied:
  - If weight is ≥1000 Da then even a small number of amide bonds or at least one sugar unit support macromolecular classification.
  - If weight is between 800 and 1000 Da, we require at least 4 amide bonds or ≥2 sugar units plus a sufficiently high atom count.
  - If weight is between 500 and 800 Da, we require at least 5 amide bonds or ≥3 sugar units and enough atoms.
  
Note: This heuristic still has limitations since “macromolecule” is an ambiguous concept.
  
Usage: Call is_macromolecule(smiles) and it returns a boolean plus a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    
    Heuristics used:
      1. If molecular weight ≥ 1000 Da then the molecule is considered macromolecular if it
         shows any evidence of being built from repeat units (e.g. one or more amide bonds or sugar units).
      2. For molecular weight between 800 and 1000 Da, at least 4 amide bonds or 2 sugar repeats 
         and at least 40 atoms are required.
      3. For molecular weight between 500 and 800 Da, at least 5 amide bonds or 3 sugar repeats 
         and at least 50 atoms are required.
      4. Otherwise, it is not classified as a macromolecule.
         
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is classified as a macromolecule, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Count total number of atoms.
    num_atoms = mol.GetNumAtoms()
    # Count amide bonds using the SMARTS pattern for a C(=O)N group.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    
    # Also count sugar-like repeating units.
    # We choose a common hexose pattern: it roughly matches a pyranose ring with hydroxyls.
    sugar_pattern = Chem.MolFromSmarts("O[C@H]1OC(O)C(O)C(O)C1O")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    num_sugar = len(sugar_matches)
    
    # Now apply our heuristic conditions.
    if mol_wt >= 1000:
        # If really high weight, even minimal signs of repetition qualify.
        if num_amide >= 1 or num_sugar >= 1 or num_atoms >= 40:
            reason = (f"Molecular weight is {mol_wt:.1f} Da which is high for a macromolecule. "
                      f"Found {num_amide} amide bonds and {num_sugar} sugar units, with {num_atoms} atoms overall.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da but no clear repetitive bonds were detected "
                      f"(0 amide bonds and 0 sugar units) in a {num_atoms}-atom structure.")
            return False, reason
    elif 800 <= mol_wt < 1000:
        if (num_amide >= 4 or num_sugar >= 2) and num_atoms >= 40:
            reason = (f"Molecular weight is {mol_wt:.1f} Da, and the molecule has {num_atoms} atoms "
                      f"with {num_amide} amide bonds and {num_sugar} sugar units, suggesting repeated subunits.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da but only {num_amide} amide bonds and "
                      f"{num_sugar} sugar units were detected in a {num_atoms}-atom structure; criteria for repeated units are not met.")
            return False, reason
    elif 500 <= mol_wt < 800:
        if (num_amide >= 5 or num_sugar >= 3) and num_atoms >= 50:
            reason = (f"Molecular weight is {mol_wt:.1f} Da, and the molecule has {num_atoms} atoms "
                      f"with {num_amide} amide bonds and {num_sugar} sugar units, suggesting multiple repeat units.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da with {num_atoms} atoms and only {num_amide} amide bonds "
                      f"and {num_sugar} sugar units; criteria for a macromolecule are not met.")
            return False, reason
    else:
        reason = (f"Molecular weight is {mol_wt:.1f} Da with only {num_atoms} atoms, {num_amide} amide bonds, "
                  f"and {num_sugar} sugar units; criteria for a macromolecule (high MW and repeated units) are not met.")
        return False, reason

# Example usage:
if __name__ == '__main__':
    # Test with an example peptide SMILES.
    test_smiles = "CC[C@H](C)[C@H](NC(=O)CN)C(=O)NCC(=O)N[C@@H](C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H]([C@@H](C)O)C(=O)NCC(=O)N[C@@H](CC(C)C)C(=O)N1CCC[C@H]1C(=O)N[C@@H](C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H]([C@@H](C)CC)C(=O)N"
    flag, explanation = is_macromolecule(test_smiles)
    print(flag)
    print(explanation)