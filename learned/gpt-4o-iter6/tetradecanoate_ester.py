"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is formed by the esterification of tetradecanoic acid (C13H27COOH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for tetradecanoate ester: C13H27C(=O)O
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O")  # C13 backbone followed by ester group
    
    # Check for presence of tetradecanoate ester functionality
    if not mol.HasSubstructMatch(tetradecanoate_pattern):
        return False, "No tetradecanoate ester group found"
    
    # Count occurrences of the pattern to ensure the presence of at least one tetradecanoate ester group
    ester_matches = mol.GetSubstructMatches(tetradecanoate_pattern)
    if len(ester_matches) < 1:
        return False, "Missing tetradecanoate esters, none detected"
    
    # Ensure accurate counting by exact number of carbons tied to ester functional group
    for ester_match in ester_matches:
        # Validate exactly 14 carbons (13 from the chain + 1 from ester) exist in each match
        carbon_count = sum(1 for atom_idx in ester_match if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6)
        if carbon_count != 14:
            return False, "Incorrect carbon count for tetradecanoate ester, expected 14 carbons including ester"

    return True, "Contains tetradecanoate ester functionality"