"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone polyphenol
Definition: A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.
This implementation uses characteristic flavonoid/isoflavone cores as a proxy to identify acrovestone-class compounds.
"""

from rdkit import Chem

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone-class polyphenol based on its SMILES string.
    Criteria (as a proxy for the acrovestone class):
      1. The molecule must form properly from the given SMILES.
      2. It should contain a flavonoid/isoflavone core typical to acrovestone-class polyphenols.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as acrovestone, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for characteristic cores.
    # Flavone core: represents a benzopyran-4-one system.
    flavone_smarts = "c1ccc2oc(=O)cc(c2)c1"
    # Isoflavone core: represents an isoflavone arrangement.
    isoflavone_smarts = "c1ccc(-c2oc(=O)c3ccccc3)cc1"

    flavone_pattern = Chem.MolFromSmarts(flavone_smarts)
    isoflavone_pattern = Chem.MolFromSmarts(isoflavone_smarts)
    
    # Verify that our SMARTS patterns compiled correctly.
    if flavone_pattern is None:
        return False, "Failed to compile flavone SMARTS pattern"
    if isoflavone_pattern is None:
        return False, "Failed to compile isoflavone SMARTS pattern"
    
    # Check for the presence of either the flavone or isoflavone core.
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Molecule contains a flavone core typical of acrovestone-class polyphenols"
    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Molecule contains an isoflavone core typical of acrovestone-class polyphenols"
    
    return False, "Molecule does not contain a flavonoid/isoflavone core structure typical of acrovestone"

# Example usage:
if __name__ == "__main__":
    # Test with one provided example: Isogenistein 7-glucoside
    test_smiles = "O1C(C(O)C(O)C(O)C1OC=2C=C3OC=C(C(=O)C3=C(O)C2)C=4C(O)=CC=CC4)CO"
    classification, explanation = is_acrovestone(test_smiles)
    print(f"Is acrovestone: {classification}\nReason: {explanation}")