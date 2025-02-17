"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone polyphenol.
Definition: A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.
This implementation first removes glycoside moieties using the FragmentParent function (from rdkit.Chem.MolStandardize)
and then checks for either a flavone (2-phenylchromen-4-one) or isoflavone (3-phenylchromen-4-one) core.
"""
from rdkit import Chem
# Import rdMolStandardize from the correct module
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:
    # In case the import fails, set to None so that the function can report an error
    rdMolStandardize = None

def is_acrovestone(smiles: str):
    """
    Determines if a molecule qualifies as an acrovestone-class polyphenol based on its SMILES string.
    Given that many acrovestone examples are glycosylated, this function first removes common sugar fragments
    (by taking the largest fragment via FragmentParent) and then checks for the presence of a flavone or isoflavone core.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as acrovestone, False otherwise.
        str: Explanation for the classification result.
    """
    if rdMolStandardize is None:
        return None, "rdMolStandardize module not available from rdkit.Chem.MolStandardize"
    
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove sugar moieties by obtaining the largest non-sugar fragment.
    try:
        parent = rdMolStandardize.FragmentParent(mol)
    except Exception as e:
        return False, f"Error during sugar removal: {str(e)}"
    if parent is None:
        return False, "Failed to determine the parent fragment after sugar removal"
    
    # Define SMARTS patterns for core structures.
    # Flavone core (2-phenylchromen-4-one)
    flavone_smarts = "c1ccc2oc(=O)cc(c2)c1"
    # Isoflavone core (3-phenylchromen-4-one)
    isoflavone_smarts = "O=C1C=C(c2ccccc2)Oc2ccccc12"
    
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    isoflavone_pattern = Chem.MolFromSmarts(isoflavone_smarts)
    if flavone_pattern is None or isoflavone_pattern is None:
        return False, "Failed to compile core SMARTS pattern(s)"
    
    # Check for the presence of the flavone core in the deglycosylated parent molecule.
    if parent.HasSubstructMatch(flavone_pattern):
        return True, "Deglycosylated molecule contains a flavone core typical of acrovestone-class polyphenols"
    
    # Check for the presence of the isoflavone core.
    if parent.HasSubstructMatch(isoflavone_pattern):
        return True, "Deglycosylated molecule contains an isoflavone core typical of acrovestone-class polyphenols"
    
    # If neither core is found, classify the molecule as non-acrovestone.
    return False, "Deglycosylated molecule does not contain a flavone or isoflavone core typical of acrovestone-class polyphenols"

# Example usage for testing:
if __name__ == "__main__":
    # Testing with one of the provided example SMILES (Pseudobaptigenin 7-O-laminaribioside)
    test_smiles = "O([C@H]1[C@H](O)C(O[C@@H](OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=5OCOC5C=C4)C1O)CO)[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO"
    classification, explanation = is_acrovestone(test_smiles)
    print("Is acrovestone:", classification)
    print("Reason:", explanation)