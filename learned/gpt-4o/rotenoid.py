"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    Rotenoids are defined by the presence of a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

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
    
    # Revised SMARTS pattern for tetrahydrochromeno[3,4-b]chromene skeleton
    # A more general pattern for rotenoid structure
    rotenoid_pattern = Chem.MolFromSmarts("C1OC(C2=CC3=C(C2O1)C=CC=C3)=O")
    
    # Check if the molecule has the rotenoid pattern
    if mol.HasSubstructMatch(rotenoid_pattern):
        return True, "Matches the rotenoid pattern"
    else:
        return False, "Does not match the rotenoid pattern"