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
    
    # A more detailed SMARTS pattern for the hopane skeleton based on known structures.
    # These SMILES are based on structures with a pentacyclic ring system.
    hopane_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC5C4(C3(C2)C1)CC5")

    if hopane_pattern is None:
        return False, "Error in SMARTS pattern"

    if mol.HasSubstructMatch(hopane_pattern):
        return True, "Contains hopane skeleton"
    else:
        return False, "No hopane skeleton recognized"