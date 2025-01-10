"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is any carboxylic ester with a carboxylic acid component
    derived from butyric acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # BUTYRATE ESTER Signature
    # Look for ester linkage pattern with the appropriate connectivity for butyric acid (CCCC(=O)O part of an ester)
    # This searches for general ester "R1-O-C(=O)R2" and checks the R1 part as the butyric acid
    butyrate_ester_pattern = Chem.MolFromSmarts("C(CCC)(=O)O")
    
    # Look for the butyrate ester substructure in the molecule
    if mol.HasSubstructMatch(butyrate_ester_pattern):
        return True, "Contains butyrate ester substructure"
    else:
        return False, "Does not contain butyrate ester substructure"

# Example usage:
# result, reason = is_butyrate_ester("CCC(=O)OCC(C)C")  # Example: isobutyl butyrate SMILES
# print(result, reason)