"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is defined as a compound containing exactly two hydroxy groups,
    but taking into account molecular context to avoid false positives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a simple alcohol diol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroxy group pattern (-OH)
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    # Count the number of hydroxy groups
    num_hydroxy_groups = len(hydroxy_matches)
    
    # A diol should have exactly two -OH groups
    if num_hydroxy_groups != 2:
        return False, f"Contains {num_hydroxy_groups} hydroxy groups, need exactly 2"

    # Additional check to avoid false positives
    # Check if the hydroxy groups are on two separate carbon atoms each bonded to no more than one oxygen
    hydroxy_carbon_matches = [match[0] for match in hydroxy_matches]
    carbon_match_criteria = all(
        mol.GetAtomWithIdx(carbon_idx).GetDegree() <= 3
        and len([nb for nb in mol.GetAtomWithIdx(carbon_idx).GetNeighbors() if nb.GetSymbol() == 'O']) == 1
        for carbon_idx in hydroxy_carbon_matches
    )
    if not carbon_match_criteria:
        return False, "Hydroxyl groups are part of complex groups, not simple diol"

    return True, "Contains exactly two hydroxy groups forming a diol"