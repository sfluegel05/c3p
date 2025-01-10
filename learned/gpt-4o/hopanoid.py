"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is characterized by a hopane skeleton, which consists of five interconnected rings forming a specific triterpenoid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for the hopane skeleton.
    # This pattern accounts for the interconnected rings and common stereochemistry patterns found with hopanoids.
    hopane_pattern = Chem.MolFromSmarts("C12CC3CC4CCC(C4)C3CC1CCC2")

    if hopane_pattern is None:
        return False, "Error in SMARTS pattern"

    # Check for substructure match in the molecule.
    if mol.HasSubstructMatch(hopane_pattern):
        return True, "Contains hopane skeleton"
    else:
        return False, "No hopane skeleton recognized"