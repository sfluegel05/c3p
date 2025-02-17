"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside

Definition: A beta-D-galactoside is any D-galactoside having beta-configuration
at its anomeric centre. In many examples (e.g. methyl beta-D-galactoside),
the glycosidic linkage is drawn as:
  CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O
In this improved approach we require that the fragment is a complete pyranose ring
(i.e. exactly one ring oxygen, five carbons, and the CH2OH substituent)
with the anomeric carbon drawn with beta stereochemistry ([C@@H]).
This more strict SMARTS pattern is intended to reduce false positive matches.
"""

from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    This function uses a SMARTS pattern that captures a complete beta-D-galactopyranoside 
    moiety. The pattern is designed to require:
      - a six-membered sugar ring (pyranose) with exactly one oxygen atom,
      - the anomeric center specified with beta stereochemistry ([C@@H]),
      - and the primary alcohol (CH2OH) substituent.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule appears to be a beta-D-galactoside, False otherwise.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS pattern for a beta-D-galactopyranoside moiety.
    # This pattern matches a six-membered ring with one ring oxygen and five carbons,
    # where the anomeric carbon (the one connected to an external O) is drawn as [C@@H]
    # (indicating beta linkage) and the ring contains the primary CH2OH substituent.
    #
    # The pattern below was built based on the common example for methyl beta-D-galactoside:
    # "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    #
    # Explanation of the pattern:
    #  - "O[C@@H]1" : a glycosidic oxygen bonded to an anomeric carbon with beta configuration.
    #  - "O[C@H](CO)" : the next ring atom is an oxygen-bearing carbon with a CH2OH group.
    #  - Followed by three consecutive "[C@H](O)" atoms;
    #  - "[C@H]1O" closes the ring with a hydroxyl substituent.
    #
    # Note: While some valid variations exist, this pattern is more selective than the previous ones.
    smarts_beta_gal = "O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    patt = Chem.MolFromSmarts(smarts_beta_gal)
    if patt is None:
        return None, "Error creating SMARTS pattern"
    
    # Search for the SMARTS pattern in the molecule.
    if mol.HasSubstructMatch(patt):
        return True, "Contains beta-D-galactoside moiety based on specific SMARTS matching"
    else:
        return False, "Beta-D-galactoside moiety not found"
        
# Example usage:
if __name__ == "__main__":
    # Example: methyl beta-D-galactoside (should be classified as True)
    test_smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_beta_D_galactoside(test_smiles)
    print(test_smiles, "->", result, "|", reason)