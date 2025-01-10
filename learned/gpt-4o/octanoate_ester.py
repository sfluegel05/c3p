"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester contains the ester linkage where the carboxylic acid component is octanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define octanoic acid ester pattern: O=C-O and a long carbon chain (7 carbon atoms connected linearly)
    # This will match ester groups with an 8-carbon chain, allowing more flexibility in esterifying groups
    octanoate_ester_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCC")
    
    # Alternative definition for octanoate ester linkage to ensure better pattern:
    detailed_pattern = Chem.MolFromSmarts("*C(=O)OCCCCCCC")
    
    # Check if the molecule contains the octanoic ester group
    if mol.HasSubstructMatch(detailed_pattern) or mol.HasSubstructMatch(octanoate_ester_pattern):
        return True, "Contains octanoic acid ester linkage"
    
    return False, "Does not contain octanoic acid ester linkage"