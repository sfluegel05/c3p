"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside

Definition: A beta-D-galactoside is any D-galactoside having beta-configuration at its anomeric centre.
Here we adopt a heuristic approach by looking for a pyranose ring fragment that is characteristic of
a beta-D-galactopyranoside. Many beta-D-galactosides have a glycosidic linkage depicted as an O–C( ) bond
joining the aglycone and a six-membered sugar ring whose stereochemistry corresponds to that of D-galactose.
We use two SMARTS patterns that capture common portrayals of this motif.
"""

from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    Uses a couple of SMARTS fragments that capture a beta-D-galactopyranoside moiety.
    (A beta glycoside is one in which the substituent at the anomeric carbon is on the same side as the 
    CH2OH group of the sugar ring; the SMARTS here are chosen from common examples and may not catch all cases.)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule appears to be a beta-D-galactoside, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns that are intended to capture the beta-D-galactoside sugar ring.
    # Many beta-D-galactosides (e.g. methyl beta-D-galactoside) feature a glycosidic oxygen 
    # attached to a six-membered ring with a pattern such as one of the following.
    #
    # Pattern 1: uses the chiral specification [C@@H] at the anomeric carbon.
    # Example match: "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    #
    # Pattern 2: a variant allowing the anomeric carbon to be specified as [C@H] (depending on how the sugar is
    # rendered in the SMILES). These patterns are heuristic.
    smarts_pattern1 = "O[C@@H]1O[C@H]([*])[C@H](O)[C@H](O)[C@H]1O"
    smarts_pattern2 = "O[C@H]1O[C@H]([*])[C@H](O)[C@H](O)[C@H]1O"
    
    patt1 = Chem.MolFromSmarts(smarts_pattern1)
    patt2 = Chem.MolFromSmarts(smarts_pattern2)
    
    # Check if either of the patterns is present in the molecule.
    if mol.HasSubstructMatch(patt1) or mol.HasSubstructMatch(patt2):
        return True, "Contains beta-D-galactoside moiety based on SMARTS matching"
    else:
        return False, "Beta-D-galactoside moiety not found"
        
# Example usage (uncomment to test):
# test_smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"  # methyl beta-D-galactoside
# print(is_beta_D_galactoside(test_smiles))