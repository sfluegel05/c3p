"""
Classifies: CHEBI:46722 carbonate ester
"""
from rdkit import Chem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    A carbonate ester has the general structure R-O-C(=O)-O-R', where R and R' are organyl groups, 
    or R-O-C(=O)-O-H for a hydrogen carbonate. This also covers cyclic carbonates.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the core carbonate ester SMARTS pattern: Any non-H atom - O - C(=O) - O - any non-H atom
    carbonate_pattern = Chem.MolFromSmarts("[!H][O][CX3](=[OX1])[O][!H]")
    
    # Check if the molecule has the carbonate group
    if not mol.HasSubstructMatch(carbonate_pattern):
        return False, "No carbonate group found"
    
    return True, "Contains a carbonate group with organyl groups attached."