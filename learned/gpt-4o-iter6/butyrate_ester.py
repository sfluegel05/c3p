"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is any carboxylic ester where the carboxylic acid component is derived from butyric acid.

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
    # Look for butyrate ester part: "CCCC(=O)O" (butyric acid with an ester linkage)
    # Allow for direct carbon bonding without constraint on stereo centers or additional substituents

    butyrate_ester_pattern = Chem.MolFromSmarts("CCCC(=O)O")  # Define signature for butyrate ester
    
    # If molecule has substructure matching butyrate ester pattern
    if mol.HasSubstructMatch(butyrate_ester_pattern):
        return True, "Contains butyrate ester substructure"
    else:
        return False, "Does not contain appropriate butyrate ester substructure"

# Example test
# result, reason = is_butyrate_ester("CCC(COC(CCC)=O)C")  # Example: 2-methylbutyl butanoate SMILES
# print(result, reason)