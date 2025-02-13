"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid moiety ('OC(=O)CCC') as part of an ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for butyrate ester group
    butyrate_pattern = Chem.MolFromSmarts("OC(=O)CCC")
    
    # Check if molecule has the butyrate ester substructure
    if mol.HasSubstructMatch(butyrate_pattern):
        return True, "Contains butyrate ester substructure"
    else:
        return False, "Does not contain butyrate ester substructure"

# Example usage:
# result, reason = is_butyrate_ester("O=C(OCCCCC)CCC")  # Example, replace with actual SMILES
# print(result, reason)