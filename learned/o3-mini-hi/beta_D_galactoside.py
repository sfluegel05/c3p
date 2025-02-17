"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside

Definition: A beta-D-galactoside is any D-galactoside having beta-configuration
at its anomeric centre. Many valid examples (e.g. methyl beta-D-galactoside)
display a pyranose ring where the anomeric carbon is drawn as [C@@H] (indicating beta linkage).
Moreover, D-galactopyranosides have the characteristic feature of a six-membered ring 
(with exactly one ring oxygen and five carbons), along with a CH2OH substituent attached
at the proper (typically C5) position and a distinctive inversion at C4 compared to glucose.
Because an exact definition in SMARTS is challenging, this improved approach requires that
a beta-D-galactoside fragment be present as a complete pyranose ring having the expected substituents.
Note: This approach will miss examples if the SMILES omit stereochemistry or use alternative conventions.
"""

from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    The method first parses the SMILES and then searches for a substructure
    that matches a refined SMARTS pattern. This pattern in effect requires that:
      - There is an external substituent binding (denoted by a wildcard [*])
        to an oxygen that is linked to a six-membered ring.
      - Within that ring the first (anomeric) carbon is drawn with beta-configuration ([C@@H]),
        and the ring is arranged as expected for galactopyranosides. In our pattern, 
        the ring is assumed to include one oxygen and five carbons, four of which bear –OH, 
        while one of the carbons bears a CH2OH substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule appears to contain a beta-D-galactoside fragment.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure that stereochemistry is assigned.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Improved SMARTS pattern:
    # Breakdown of the pattern:
    #   [*]               : some substituent attached to...
    #   O                 : an oxygen atom (the glycosidic oxygen)
    #   [C@@H]1         : an anomeric carbon with beta stereochemistry initiating ring '1'
    #   [C@H](O)        : a ring carbon bearing a hydroxyl group
    #   [C@H](O)        : a second ring carbon bearing a hydroxyl group
    #   [C@@H](O)       : a third ring carbon bearing a hydroxyl group (configuration inverted as for galactose)
    #   [C@H](CO)       : a fourth ring carbon bearing the CH2OH substituent
    #   O1              : a ring oxygen that closes the ring.
    #
    # This pattern is selective for a complete beta-D-galactopyranoside-fragment.
    smarts_beta_gal = "[*]O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](CO)O1"
    patt = Chem.MolFromSmarts(smarts_beta_gal)
    if patt is None:
        return None, "Error creating SMARTS pattern"
    
    # Look for the substructure match in the molecule.
    if mol.HasSubstructMatch(patt):
        return True, "Contains beta-D-galactoside moiety based on refined SMARTS matching"
    else:
        return False, "Beta-D-galactoside moiety not found"

# Example usage:
if __name__ == "__main__":
    # Test with methyl beta-D-galactoside (should be classified as True)
    test_smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_beta_D_galactoside(test_smiles)
    print(test_smiles, "->", result, "|", reason)