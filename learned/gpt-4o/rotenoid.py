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
    
    # Revised SMARTS pattern for rotenoid core structure
    # A simplified and possibly inclusive pattern for the tetrahydrochromeno[3,4-b]chromene core
    rotenoid_core_pattern = Chem.MolFromSmarts("C1Oc2ccccc2Oc3ccccc13")
    
    # Check if the molecule contains the reduced rotenoid core pattern
    if mol.HasSubstructMatch(rotenoid_core_pattern):
        return True, "Matches the rotenoid core pattern"
    else:
        return False, "Does not match the rotenoid core pattern"