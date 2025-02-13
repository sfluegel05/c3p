"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids are characterized by a C6-C3-C6 carbon skeleton attached to a heterocyclic benzopyran ring.
    They include subclasses such as flavones, flavonols, isoflavones, and flavanones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for flavonoid subclasses
    flavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)cc2")  # Basic flavone structure
    
    # Other flavonoid patterns can be added here for more comprehensive detection...
    flavonol_pattern = Chem.MolFromSmarts("c1ccc2c(c1)c(=O)c3cc(O)ccc3o2")
    isoflavone_pattern = Chem.MolFromSmarts("c1ccc2c(c1)cc(=O)oc2c3ccc(c(c3)O)O")

    # Check for matching any of the flavonoid patterns
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Identified as flavone structure"
    
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Identified as flavonol structure"
    
    if mol.HasSubstructMatch(isoflavone_pattern):
        return True, "Identified as isoflavone structure"

    # Could add more checks for other subclasses of flavonoids

    return False, "No flavonoid pattern matched"