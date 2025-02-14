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

    # Define enriched SMARTS patterns for identifying carbamate esters
    carbamate_ester_patterns = [
        Chem.MolFromSmarts("COC(=O)N"),          # Simple carbamate ester
        Chem.MolFromSmarts("COC(=O)N(C)"),       # N-methyl carbamate
        Chem.MolFromSmarts("COC(=O)N(C)C"),      # N,N-dimethyl carbamate
        Chem.MolFromSmarts("N(C)(C)C(=O)O"),     # Tertiary amine carbamate
        Chem.MolFromSmarts("O=C(O)N([#6])"),     # N-substituted carbamate with ester linkage
        Chem.MolFromSmarts("O=C([O,N])[N]"),     # General pattern for carbamate linkage
    ]

    # Validate carbamate ester presence using structured substructure searches
    for pattern in carbamate_ester_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carbamate ester functional group"

    return False, "Does not contain carbamate ester functional group"