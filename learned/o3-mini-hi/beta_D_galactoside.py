"""
Classifies: CHEBI:28034 beta-D-galactoside
"""
"""
Classifies: beta-D-galactoside

Definition: Any D-galactoside having beta-configuration at its anomeric centre.
A beta-D-galactoside is defined here as a molecule that contains a complete beta-D-galactopyranoside
fragment. In our approach we require that an external substituent is linked via an oxygen atom 
to the anomeric carbon (which must be in beta-configuration, represented as [C@@H]), and that 
this anomeric carbon is part of a six-membered ring consisting of one ring oxygen and five carbons,
with the expected substituents (three hydroxyl groups and one CH2OH group at the proper position).
Because galactose differs from other hexopyranoses in its stereochemistry (typically an inversion
at C3 or C4 compared to glucose), our SMARTS is designed to enforce the expected order.
Note: This approach is strict and may fail if stereochemistry is omitted or if alternative 
representations are used.
"""

from rdkit import Chem

def is_beta_D_galactoside(smiles: str):
    """
    Determines if a molecule is a beta-D-galactoside based on its SMILES string.
    
    This function first parses the SMILES string and then searches for a specific substructure 
    that represents a beta-D-galactopyranoside fragment. The pattern requires that:
      - An external substituent is attached via an oxygen ([*]O)
      - That oxygen links to an anomeric carbon drawn with beta configuration ([C@@H])
      - The anomeric carbon is part of a six-membered ring arranged as:
            O (ring oxygen) - [C@H](CO) - [C@H](O) - [C@H](O) - [C@H] (ring closure with –OH)
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule contains a beta-D-galactoside fragment, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Force RDKit to assign stereochemistry.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Our refined SMARTS now requires a complete beta-D–galactopyranoside fragment:
    #   [*]O -> any substituent linked via an oxygen to...
    #   [C@@H]1 -> the anomeric carbon in beta configuration (starting a ring '1')
    #   O -> ring oxygen (the next atom in the pyranose ring)
    #   [C@H](CO) -> first ring carbon bearing a CH2OH substituent (typical for C2 in D-galactose)
    #   [C@H](O)  -> second ring carbon bearing an -OH
    #   [C@H](O)  -> third ring carbon bearing an -OH
    #   [C@H]1O  -> fourth ring carbon closing the ring with an -OH substituent.
    #
    # This fragment exactly represents a six-membered pyranose where the linkage formation is explicit.
    smarts_beta_gal = "[*]O[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    patt = Chem.MolFromSmarts(smarts_beta_gal)
    if patt is None:
        return None, "Error creating SMARTS pattern"
    
    # Look for a substructure match with the refined and strict pattern.
    if mol.HasSubstructMatch(patt):
        return True, "Contains beta-D-galactoside moiety based on refined SMARTS matching"
    else:
        return False, "Beta-D-galactoside moiety not found"

# Example usage (for testing):
if __name__ == "__main__":
    # Testing with a confirmed beta-D-galactoside: methyl beta-D-galactoside.
    test_smiles = "CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O"
    result, reason = is_beta_D_galactoside(test_smiles)
    print(test_smiles, "->", result, "|", reason)