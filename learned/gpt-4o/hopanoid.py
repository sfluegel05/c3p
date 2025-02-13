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
    
    # A detailed SMARTS pattern for the hopane skeleton.
    # We are attempting to define rings and sterochemistry known in hopanoids
    hopane_pattern = Chem.MolFromSmarts("C1CCC2C3C4[C@H]5CC[C@@H](C5)CC4CCC3CC2CCC1")

    if hopane_pattern is None:
        return False, "Error in SMARTS pattern"

    # Check for substructure match in the molecule
    if mol.HasSubstructMatch(hopane_pattern):
        return True, "Contains hopane skeleton"
    else:
        return False, "No hopane skeleton recognized"