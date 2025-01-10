"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is characterized by a five-membered lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced SMARTS pattern for gamma-lactone (five-membered lactone ring)
    # This pattern allows for flexible substituents on the lactone ring
    gamma_lactone_pattern = Chem.MolFromSmarts("[#6]-1-[O;$(C=O)]-C(=O)-[C;R]1")
    
    # Check for gamma-lactone pattern
    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains a five-membered lactone ring"

    return False, "Does not contain a five-membered lactone ring"