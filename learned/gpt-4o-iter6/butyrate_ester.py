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

    # The SMARTS pattern for a butyrate ester aims at capturing structures like:
    # CCCC(=O)O , where the exact ester structure is derived from butyric acid.
    # We use "*" to denote variably connected atoms.

    # Butyric acid part: "O=C(O)CCC" (ester has replaced the carboxylic O-H)
    butyrate_ester_pattern = Chem.MolFromSmarts("*O=C(OCCC)")

    # Check if the molecule has this substructure pattern
    if mol.HasSubstructMatch(butyrate_ester_pattern):
        return True, "Contains a butyrate ester substructure"
    else:
        return False, "Does not contain appropriate butyrate ester substructure"
    
# Example test
# result, reason = is_butyrate_ester("CCC(COC(CCC)=O)C")  # Example: 2-methylbutyl butanoate SMILES
# print(result, reason)