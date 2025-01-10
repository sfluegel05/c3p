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
    
    # Define a more flexible SMARTS pattern for butyrate ester group
    # Look for -C-C(=O)O-C structure where the C before O could vary
    butyrate_pattern = Chem.MolFromSmarts("CCC(=O)O")
    
    # Ensure the butanoate is part of a larger ester linkage
    esters_side_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    
    # Combine general butyrate pattern with ester linkage check
    if mol.HasSubstructMatch(butyrate_pattern) and mol.HasSubstructMatch(esters_side_pattern):
        return True, "Contains butyrate ester substructure"
    else:
        return False, "Does not contain butyrate ester substructure"

# Example usage:
# result, reason = is_butyrate_ester("CCC(=O)OC")  # Example: isobutyl butyrate SMILES
# print(result, reason)