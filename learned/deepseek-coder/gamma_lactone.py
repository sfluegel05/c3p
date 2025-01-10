"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: CHEBI:37671 gamma-lactone
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a five-membered lactone ring.

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

    # Define a more flexible gamma-lactone pattern: a five-membered ring with one oxygen and a carbonyl group
    gamma_lactone_pattern = Chem.MolFromSmarts("[O;R1]1[C;R1](=O)[C;R1][C;R1][C;R1]1")
    
    # Check if the molecule matches the gamma-lactone pattern
    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains a five-membered lactone ring (gamma-lactone)"
    else:
        return False, "No five-membered lactone ring found"