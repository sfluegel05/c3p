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

    # Define the flavone core using SMARTS pattern for 2-phenylchromen-4-one
    flavone_core_pattern = Chem.MolFromSmarts('O=C1C=CC(=CC1=O)c1ccccc1')

    # Check if the molecule contains the flavone core structure
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No flavone core structure found"

    return True, "Contains flavone core structure (2-phenylchromen-4-one)"

# Example usage:
# result, reason = is_flavones("COc1cc(O)c2c(c1)oc(cc2=O)-c1ccc(O)cc1")
# print(result, reason)