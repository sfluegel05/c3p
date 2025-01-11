"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester contains an ester linkage with an octanoic acid part (8-carbon chain).
    
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

    # Look for ester bond pattern (R-C(=O)O-R')
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester bond found"

    # Check for 8-carbon chain linked to ester group
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    if not mol.HasSubstructMatch(octanoate_pattern):
        return False, "No octanoate chain found associated with ester bond"

    return True, "Molecule contains an octanoate ester group"

# Example Call
# result, reason = is_octanoate_ester("O(C(CCCCCCC)=O)C=1C=C2C=CC=CC2=CC1")  # Example SMILES for testing
# print(result, reason)