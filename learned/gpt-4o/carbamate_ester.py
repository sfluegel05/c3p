"""
Classifies: CHEBI:23003 carbamate ester
"""
from rdkit import Chem

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as any ester of carbamic acid or its N-substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbamate ester, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for typical carbamate ester forms
    # This pattern might need to match various typical bond configurations.
    carbamate_ester_patterns = [
        Chem.MolFromSmarts("COC(=O)N"),  # Methyl ester variant
        Chem.MolFromSmarts("COC(=O)N(C)"),  # Methyl N-alkyl carbamate
        Chem.MolFromSmarts("COC(=O)NC"),    # N-substituted, possible secondary
        Chem.MolFromSmarts("CC(=O)OC(=O)N"),  # E.g., ethyl/other esters
    ]
    
    # Check if the molecule has any of the carbamate ester substructures
    for pattern in carbamate_ester_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbamate ester functional group"

    return False, "Does not contain carbamate ester functional group"

# Note: Pass a SMILES string to the is_carbamate_ester() function to classify.