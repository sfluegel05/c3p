"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone polyphenol.
Definition: A polyphenol from Acronychia pedunculata with moderate antioxidant and antityrosinase activities.
This implementation attempts to deglycosylate the molecule using rdMolStandardize.FragmentParent and then check for
a flavone, isoflavone, or generic chromen-4-one core.
"""
from rdkit import Chem
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:
    rdMolStandardize = None

def is_acrovestone(smiles: str):
    """
    Determines if a molecule qualifies as an acrovestone-class polyphenol based on its SMILES string.
    Since many acrovestone examples are glycosylated, this function first tries to remove sugar moieties
    (via FragmentParent) and then checks if the remaining core contains a flavone (2-phenylchromen-4-one),
    an isoflavone (3-phenylchromen-4-one), or a more generic chromen-4-one structure.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as acrovestone, False otherwise.
        str: Explanation for the classification.
    """
    if rdMolStandardize is None:
        return None, "rdMolStandardize module not available from rdkit.Chem.MolStandardize"
    
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempt to remove glycoside (sugar) moieties by obtaining the parent fragment.
    try:
        parent = rdMolStandardize.FragmentParent(mol)
    except Exception as e:
        return False, f"Error during sugar removal: {str(e)}"
    if parent is None:
        return False, "Failed to determine the parent fragment after sugar removal"
    
    # Define SMARTS patterns for core structures.
    # Pattern 1: Flavone core (2-phenylchromen-4-one)
    flavone_smarts = "c1ccc2oc(=O)cc(c2)c1"
    # Pattern 2: Isoflavone core (3-phenylchromen-4-one)
    isoflavone_smarts = "O=C1C=C(c2ccccc2)Oc2ccccc12"
    # Pattern 3: Broad chromen-4-one core pattern (any fused benzopyrone ring system)
    generic_chromen_smarts = "O=C1C=CC2=CC=CC=C2O1"
    
    # Compile SMARTS patterns.
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    isoflavone_pattern = Chem.MolFromSmarts(isoflavone_smarts)
    generic_pattern = Chem.MolFromSmarts(generic_chromen_smarts)
    if flavone_pattern is None or isoflavone_pattern is None or generic_pattern is None:
        return False, "Failed to compile one or more SMARTS patterns for core detection"
    
    # Check if the deglycosylated molecule has one of the expected cores.
    if parent.HasSubstructMatch(flavone_pattern):
        return True, "Deglycosylated molecule contains a flavone core typical of acrovestone-class polyphenols"
    if parent.HasSubstructMatch(isoflavone_pattern):
        return True, "Deglycosylated molecule contains an isoflavone core typical of acrovestone-class polyphenols"
    if parent.HasSubstructMatch(generic_pattern):
        return True, "Deglycosylated molecule contains a chromen-4-one core typical of acrovestone-class polyphenols"
    
    # If none of the core patterns are found, then the molecule is not in the acrovestone class.
    return False, "Deglycosylated molecule does not contain a recognizable flavone, isoflavone, or chromen-4-one core typical of acrovestone-class polyphenols"

# Example usage for testing:
if __name__ == "__main__":
    # Testing with one of the provided example SMILES (Pseudobaptigenin 7-O-laminaribioside)
    test_smiles = "O([C@H]1[C@H](O)C(O[C@@H](OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=5OCOC5C=C4)C1O)CO)[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO"
    classification, explanation = is_acrovestone(test_smiles)
    print("Is acrovestone:", classification)
    print("Reason:", explanation)