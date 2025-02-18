"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside with a galactose head group.
A galactosylceramide is a glycosphingolipid consisting of a ceramide (sphingoid base + fatty acid 
attached via an amide bond) linked to a monosaccharide that is galactose.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    The algorithm looks for:
      1. A galactose-containing sugar ring (in beta-D configuration) as the head group.
      2. A ceramide motif (an amide bond linking a sphingoid base and a fatty acid chain).
      3. An overall molecular size consistent with glycosphingolipids.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a galactosylceramide, else False.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Step 1. Look for a galactose-based sugar moiety ---
    # A typical beta-D-galactopyranoside substructure often appears with a glycosidic bond.
    # Here we use two patterns:
    # Pattern A: Includes an initial CH2 from the glycosidic linkage
    sugar_pattern1 = Chem.MolFromSmarts("CO[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    # Pattern B: The sugar ring on its own.
    sugar_pattern2 = Chem.MolFromSmarts("O[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    
    has_sugar = False
    if sugar_pattern1 is not None and mol.HasSubstructMatch(sugar_pattern1):
        has_sugar = True
    elif sugar_pattern2 is not None and mol.HasSubstructMatch(sugar_pattern2):
        has_sugar = True
        
    if not has_sugar:
        return False, "No galactosyl sugar moiety detected."
    
    # --- Step 2. Look for a ceramide-like motif ---
    # We check for an amide bond which is present in ceramides.
    # (Note: many molecules have amide groups, so this is a simple but not perfect check.)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if amide_pattern is None or not mol.HasSubstructMatch(amide_pattern):
        return False, "No ceramide-like amide bond detected."
    
    # Optionally, one may try to enforce that the fatty acyl chain attached to the amide is long.
    # For example, check that the molecule has many heavy atoms (typical glycosphingolipids are large).
    if mol.GetNumHeavyAtoms() < 30:
        return False, "Molecule appears too small to be a galactosylceramide."
    
    # --- Additional checks could be implemented here, such as verifying the long aliphatic chain lengths,
    # ensuring a sphingoid backbone, or even using stereochemical criteria. 
    # However, the substructure queries above capture the key aspects of galactosylceramides.

    return True, "Contains a galactosyl sugar moiety and a ceramide amide bond, consistent with a galactosylceramide."

# Example usage (uncomment for testing):
# smiles_example = "CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_galactosylceramide(smiles_example)
# print(result, reason)