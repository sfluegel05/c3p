"""
Classifies: CHEBI:24043 flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    Flavones are flavonoids with a specific 2-phenylchromen-4-one core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into a RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavone core using a more flexible SMARTS pattern
    # This accommodates modifications and typical substitution positions
    flavone_core_pattern = Chem.MolFromSmarts('O=C1C=CC(=CO1)c2ccccc2')  # Main core with possible phenyl substitutions
        
    # Check if the molecule contains the flexible flavone core structure
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No flavone core structure found (even with variability)"

    return True, "Contains flavone core structure (with 2-phenylchromen-4-one)"

# We can test it with some known flavone SMILES strings
# result, reason = is_flavones("COc1ccc(cc1)-c1cc(=O)c2c(O)cc(O)cc2o1")
# print(result, reason)