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
    
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a gamma-lactone, which is a five-membered ring
    # with a lactone group (cyclic ester), specifically: [O]-C-C-C-C(=O)
    gamma_lactone_pattern = Chem.MolFromSmarts("C1COC(=O)C1")
    
    # Check for gamma-lactone pattern
    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains a five-membered gamma-lactone ring"
    
    return False, "Does not contain a five-membered gamma-lactone ring"