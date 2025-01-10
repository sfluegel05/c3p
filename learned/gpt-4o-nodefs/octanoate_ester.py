"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester contains an ester linkage where the alcohol moiety
    has formed an ester with octanoic acid, characterized by an 8-carbon aliphatic chain.

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

    # Define an ester linkage pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for basic ester bond
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond found"
  
    # Define an octanoate chain pattern via a dentistive count of the carbon chain on ester
    # Using a pattern of 8 contiguous carbon atoms as necessary for octanoate
    octanoate_chain_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O[*]")
    
    # Check if there's a match to the octanoate chain
    if mol.HasSubstructMatch(octanoate_chain_pattern):
        return True, "Molecule contains an octanoate ester group"

    # If no match found for specific octanoate ester composition  
    return False, "No octanoate ester group with 8-carbon chain found"