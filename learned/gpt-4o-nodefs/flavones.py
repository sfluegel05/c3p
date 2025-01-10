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

    # Define a flexible flavone core pattern
    flavone_core_pattern = Chem.MolFromSmarts('c1cc2oc(=O)cc(c2c1)-c1ccccc1')  # flexible flavone core
    
    # Check if the molecule contains the flexible flavone core structure
    if not mol.HasSubstructMatch(flavone_core_pattern):
        return False, "No flavone core structure found"

    return True, "Contains flavone core structure (2-phenylchromen-4-one)"

# The function can now be tested with known flavone and non-flavone SMILES strings.
# Example usage with flavone SMILES:
# result, reason = is_flavones("COc1cc(ccc1)-c1cc(=O)c2c(O)c(O)cc(O)c2o1")
# print(result, reason)