"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    A rotenoid consists of a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton and its derivatives.

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

    # Updated SMARTS pattern to represent the rotenoid core structure.
    # This pattern highlights the common polycyclic structure of rotenoids.
    rotenoid_pattern = Chem.MolFromSmarts(
        "O1C2=C(C=3C(OC2)=CC=4OCOC4=C3)C(=O)C5=C1C=C6OC=CC6=C5"
    )

    # Ensure the pattern determines core features of rotenoids.
    if rotenoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Perform the substructure match
    if not mol.HasSubstructMatch(rotenoid_pattern):
        return False, "No rotenoid core structure found"

    return True, "Contains the rotenoid core structure"