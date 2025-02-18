"""
Classifies: CHEBI:2440 acrovestone
"""
"""
Classifies: Acrovestone polyphenol
Definition: A polyphenol that is isolated from Acronychia pedunculata and exhibits moderate antioxidant and antityrosinase activities.
This implementation searches for characteristic flavonoid/isoflavone cores as a proxy for the acrovestone class.
"""

from rdkit import Chem

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone-class polyphenol based on its SMILES string.
    Criteria (as a proxy for the acrovestone class):
      1. The molecule must be valid.
      2. It should contain a polyphenolic flavonoid core (e.g. a benzopyran-4-one core) or an isoflavone core.
         Many acrovestone examples (such as various daidzein, genistein or biochanin A glycosides) contain such cores.
    
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
    
    # Define SMARTS patterns for a flavone core and an isoflavone core.
    # The flavone core (2-phenyl-4H-1-benzopyran-4-one) is represented by the following pattern.
    # This pattern should match a fused benzene and pyran-4-one ring.
    flavone_pattern = Chem.MolFromSmarts("c1ccc2oc(=O)cc(c2)c1")
    # The isoflavone core (where the B-ring is attached at position 3) is represented by this pattern.
    isoflavone_pattern = Chem.MolFromSmarts("c1ccc(-c2oc(=O)c3ccccc3)cc1")
    
    # Check if the molecule contains either the flavone or isoflavone substructure.
    if mol.HasSubstructMatch(flavone_pattern) or mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Molecule contains a flavonoid/isoflavone core typical of acrovestone-class polyphenols"
    else:
        return False, "Molecule does not contain a flavonoid/isoflavone core structure typical of acrovestone"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples:
    test_smiles = "O1C(C(O)C(O)C(O)C1OC=2C=C3OC=C(C(=O)C3=C(O)C2)C=4C(O)=CC=CC4)CO"  # Isogenistein 7-glucoside
    is_class, reason = is_acrovestone(test_smiles)
    print(f"Is acrovestone: {is_class}\nReason: {reason}")