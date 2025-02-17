"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone polyphenol
Definition: A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.
This implementation uses the presence of a characteristic flavone or isoflavone core (after removing glycoside moieties)
as a proxy for acrovestone-class compounds.
"""
from rdkit import Chem
from rdkit.Chem import rdMolStandardize

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone-class polyphenol based on its SMILES string.
    As many acrovestone examples are glycosylated derivatives, this method first removes sugar fragments
    and then searches for either a flavone (2-phenylchromen-4-one) or an isoflavone (3-phenylchromen-4-one)
    substructure.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as acrovestone, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES; return error message for invalid input
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove sugar moieties by taking the largest (non-sugar) fragment.
    # FragmentParent removes common sugar fragments.
    parent = rdMolStandardize.FragmentParent(mol)
    if parent is None:
        return False, "Could not determine parent molecule after sugar removal"
    
    # Define SMARTS patterns for cores.
    # Flavone core (2-phenylchromen-4-one)
    flavone_smarts = "c1ccc2oc(=O)cc(c2)c1"
    # Isoflavone core (3-phenylchromen-4-one)
    isoflavone_smarts = "O=C1C=C(c2ccccc2)Oc2ccccc12"
    
    # Compile SMARTS patterns.
    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    isoflavone_pattern = Chem.MolFromSmarts(isoflavone_smarts)
    if flavone_pattern is None:
        return False, "Failed to compile flavone SMARTS pattern"
    if isoflavone_pattern is None:
        return False, "Failed to compile isoflavone SMARTS pattern"
    
    # Check if the deglycosylated molecule contains a flavone or isoflavone core.
    if parent.HasSubstructMatch(flavone_pattern):
        return True, "Molecule (after sugar removal) contains a flavone core typical of acrovestone-class polyphenols"
    if parent.HasSubstructMatch(isoflavone_pattern):
        return True, "Molecule (after sugar removal) contains an isoflavone core typical of acrovestone-class polyphenols"
    
    # If neither core is found, then classify as non-acrovestone.
    return False, "Molecule does not contain a flavone or isoflavone core typical of acrovestone-class polyphenols"

# Example usage for testing:
if __name__ == "__main__":
    # Example: using one of the provided SMILES (Pseudobaptigenin 7-O-laminaribioside)
    test_smiles = "O([C@H]1[C@H](O)C(O[C@@H](OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=5OCOC5C=C4)C1O)CO)[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO"
    classification, explanation = is_acrovestone(test_smiles)
    print(f"Is acrovestone: {classification}\nReason: {explanation}")