"""
Classifies: CHEBI:26125 phytosterols
"""
#!/usr/bin/env python
"""
Classifies: Phytosterols – Sterols similar to cholesterol that occur in plants.
Phytosterols are characterized by a central steroid nucleus 
(a cyclopentanoperhydrophenanthrene core: three six–membered rings and one five–membered ring)
and often a hydroxyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule belongs to the phytosterols class based on its SMILES string.
    
    In this improved implementation we:
      1. Parse the SMILES string.
      2. Use a SMARTS pattern that attempts to capture the steroid nucleus 
         (i.e. the fused ring system of 3 six‐membered rings and 1 five–membered ring).
      3. Check for at least one hydroxyl (–OH) group.
      4. Ensure the molecular weight is within a typical sterol range (roughly 300–600 Da).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phytosterol, False otherwise.
        str: Explanation (reason) for the classification decision.
    """
    # Parse the SMILES input
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the steroid nucleus (a simplified cyclopentanoperhydrophenanthrene core)
    # This pattern looks for a typical steroid nucleus: three 6-membered rings and one 5-membered ring fused together.
    # Note: While this pattern does not capture every possible stereochemical or double-bond variation,
    # it is useful as a fast filter.
    steroid_smarts = "[$([C@@H]1CC[C@]2(C)CC[C@H]3[C@@H]([C@@H]1)CC[C@@H]23)]"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Error creating steroid nucleus SMARTS pattern"
    
    # Check if the molecule contains the steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Molecule does not contain a characteristic steroid nucleus"
    
    # Check for the presence of at least one hydroxyl group (–OH).
    # Most phytosterols have a 3β-hydroxyl group or related oxygen functionality.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxyl_pattern is None:
        return False, "Error creating hydroxyl SMARTS pattern"
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Molecule lacks a hydroxyl group typical of phytosterols"
    
    # Verify the molecular weight is within a typical sterol range (~300-600 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the typical sterol range"
    
    # If all conditions are met, classify as a phytosterol.
    return True, "Molecule contains a steroid nucleus with a hydroxyl group and an appropriate molecular weight, typical of phytosterols"

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Testing with one of the provided examples: ergosta-5,7-dien-3beta-ol
    test_smiles = "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](C)C(C)C"
    result, reason = is_phytosterols(test_smiles)
    print(result, reason)