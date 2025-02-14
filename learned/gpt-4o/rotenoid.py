"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    A rotenoid consists of a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Revised SMARTS pattern for rotenoid structure
    # The pattern aims to capture essential features of the rotenoid scaffold, including the cis-fused rings.
    rotenoid_pattern = Chem.MolFromSmarts("C12OC3=C(O1)C=CC4=C3OC5=CC=CC5=C24")
    
    # Check if the molecule has the rotenoid pattern
    if mol.HasSubstructMatch(rotenoid_pattern):
        return True, "Matches the rotenoid pattern"
    else:
        return False, "Does not match the rotenoid pattern"