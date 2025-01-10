"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined by a five-membered lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to create a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the five-membered lactone ring.
    gamma_lactone_pattern = Chem.MolFromSmarts("C1OC(=O)CC1")
    
    # Check for the presence of a gamma-lactone ring using the defined pattern.
    if mol.HasSubstructMatch(gamma_lactone_pattern):
        return True, "Contains gamma-lactone ring"
    
    return False, "No gamma-lactone ring found"

__metadata__ = {
   'chemical_class': {
       'name': 'gamma-lactone',
       'definition': 'A lactone having a five-membered lactone ring.',
   }
}