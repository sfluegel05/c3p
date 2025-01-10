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

    # Update the SMARTS pattern to better match the rotenoid core
    # This pattern attempts to capture the tetracyclic core structure of rotenoids
    # Adjusts for cis-fusion via precise sequence of aromatic and aliphatic carbons
    rotenoid_pattern = Chem.MolFromSmarts("Oc1cc2ccc3c(c2c(c1)C(=O)c4ccc5ccc(OC)cc5c43)")

    if rotenoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Match against the SMARTS pattern
    if not mol.HasSubstructMatch(rotenoid_pattern):
        return False, "No tetracyclic rotenoid core structure found"
    
    # Further checks could evaluate substitutions and stereo-chemical aspects
    # This typically requires empirical data on structure-activity relationships
    
    return True, "Contains the tetracyclic rotenoid core structure"

# Testing with known examples outside the function tests