"""
Classifies: CHEBI:32955 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is a cyclic ether with an oxygen atom in a three-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxide, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the epoxide SMARTS pattern (three-membered cyclic ether)
    epoxide_pattern = Chem.MolFromSmarts("C1OC1")

    # Check for the epoxide structure
    if mol.HasSubstructMatch(epoxide_pattern):
        return True, "Epoxide substructure found"
    else:
        return False, "No epoxide substructure found"

# Example usage
# print(is_epoxide("C1OC1"))  # Simple epoxide example