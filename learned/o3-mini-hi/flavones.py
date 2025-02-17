"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones
A flavone is a member of the flavonoid class with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton.
This function checks whether a given SMILES string has that core structure.
"""

from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone (or a substituted derivative) based on its SMILES string.
    The flavone core is defined as a 2-phenylchromen-4-one skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains the flavone skeleton, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the 2-phenylchromen-4-one (flavone) core.
    # This pattern represents:
    #  - a benzene ring (c1ccc...),
    #  - fused to a heterocycle containing an oxygen and a ketone (chromen-4-one),
    #  - with an attached phenyl group at position 2.
    # The SMARTS pattern below is:
    #   c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3
    # This should match the core of flavones and their substituted derivatives.
    flavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(c2)-c3ccccc3")
    if flavone_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Check if the molecule contains the flavone core substructure.
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains the flavone (2-phenylchromen-4-one) skeleton."
    else:
        return False, "Does not contain a 2-phenylchromen-4-one (flavone) core."

# For demonstration (remove or comment out these lines if using as an imported module):
if __name__ == "__main__":
    # Test with one of the provided examples e.g. 3-methoxyapigenin:
    test_smiles = "COc1c(oc2cc(O)cc(O)c2c1=O)-c1ccc(O)cc1"
    result, reason = is_flavones(test_smiles)
    print("Test molecule:", test_smiles)
    print("Is flavone?", result)
    print("Reason:", reason)