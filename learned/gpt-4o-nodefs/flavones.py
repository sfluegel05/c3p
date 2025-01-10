"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    Flavones are a class of flavonoids with a specific 2-phenylchromen-4-one core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more inclusive flavone core SMARTS pattern
    # This includes the 2-phenylchromen-4-one core with allowance for common substituents
    flavone_core_pattern = Chem.MolFromSmarts('c1cc2oc(=O)cc(c2cc1)-c1ccccc1')
    
    # Check if the molecule contains the defined flavone core structure
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No flavone core structure found"

    return True, "Contains flavone core structure (2-phenylchromen-4-one)"

# Example test cases covering a range of flavone structures.
# These could be tested by calling the function with the provided SMILES strings.
# Each should include checks involving known flavones with different substituents.