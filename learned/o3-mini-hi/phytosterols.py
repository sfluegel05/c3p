"""
Classifies: CHEBI:26125 phytosterols
"""
#!/usr/bin/env python
"""
Classifies: Phytosterols – Sterols similar to cholesterol that occur in plants.
Phytosterols exhibit a four‐ring (cyclopentanoperhydrophenanthrene) steroid nucleus and in most cases feature a hydroxyl group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule belongs to the phytosterols class based on its SMILES string.
    Phytosterols are plant sterols that contain the steroid nucleus
    (a fused system of three six-membered rings and one five-membered ring) and
    usually have a hydroxyl group (typically at the 3β position).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a phytosterol, False otherwise.
        str: Explanation (reason) for the classification decision.
    """

    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for identifying a steroid nucleus.
    # This approximate pattern is designed to match the tetracyclic system 
    # (three six-membered rings fused to one five-membered ring)
    steroid_nucleus_smarts = "[C@H]1CC[C@]2(CC[C@@H]1C)CCC3C2CC[C@H]3"
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if steroid_nucleus is None:
        return False, "Error creating steroid nucleus SMARTS pattern"
    
    # Check for a steroid nucleus substructure in molecule
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "Molecule does not contain the characteristic steroid nucleus"

    # Check for presence of at least one hydroxyl group (–OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxyl_pattern is None:
        return False, "Error creating hydroxyl SMARTS pattern"
    
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Molecule lacks a hydroxyl group (characteristic of sterols)"

    # Verify the molecular weight is within a typical range for sterols (~300-600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the typical range for phytosterols"

    return True, "Molecule has a steroid nucleus and a hydroxyl group typical of phytosterols"

# Example usage (if you wish to test):
if __name__ == "__main__":
    test_smiles = "OC1CC2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C"  # cholestan-3-ol
    result, reason = is_phytosterols(test_smiles)
    print(result, reason)