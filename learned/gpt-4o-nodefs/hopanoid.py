"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    Hopanoids are a class of tetracyclic triterpenoids.

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

    # Define a hopane backbone pattern for basic hopanoids
    hopane_pattern = Chem.MolFromSmarts("C1CC2[C@]3([C@]4([C@]5(C(CC4)C(CC5)C)CC3)C)C2CC1")
    if mol.HasSubstructMatch(hopane_pattern):
        return True, "Contains hopane-like backbone"

    return False, "No hopane-like backbone found"

# Examples given should match this pattern or a close variation with functional groups