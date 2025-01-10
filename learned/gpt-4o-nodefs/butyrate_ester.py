"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester contains a butyric acid moiety (4-carbon chain) attached via an ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the butyrate ester pattern: R'C(=O)O-R, with specific pattern for butyrate (4-carbon chain)
    butyrate_ester_pattern = Chem.MolFromSmarts("C(=O)OCC[CH2]C")
    
    # Check for both main and branched variations of the butyrate plus ester linkage
    standard_matches = mol.GetSubstructMatches(butyrate_ester_pattern)
    
    # Also check for common branched or isomer variations of butyrate esters
    branched_butyrate_pattern = Chem.MolFromSmarts("C(=O)O[C](C)CC")
    branched_matches = mol.GetSubstructMatches(branched_butyrate_pattern)

    if standard_matches or branched_matches:
        return True, "Contains butyrate ester functional group with valid structure"

    return False, "No valid butyrate ester group found"