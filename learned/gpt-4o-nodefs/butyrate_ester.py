"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    Includes more pattern variations and structure verification to capture true butyrate esters.

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

    # Define SMARTS patterns for butyrate ester variations, considering both straight-chain and branched forms
    ester_linkage_pattern = "C(=O)O"
    carbon_chain_4 = "CCCC"
    branched_chain = "C(C)CC"
    ester_patterns = [
        Chem.MolFromSmarts(ester_linkage_pattern + carbon_chain_4),  # Straight-chain butyrate
        Chem.MolFromSmarts(ester_linkage_pattern + branched_chain),  # Branched-chain butyrate
        Chem.MolFromSmarts(ester_linkage_pattern + "C[CH2]C"),       # Account for middle carbon variation
    ]
    
    # Check for matches with any of the ester pattern variations
    for pattern in ester_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains butyrate ester functional group with valid structure"
    
    return False, "No valid butyrate ester group found"