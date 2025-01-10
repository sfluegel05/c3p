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

    # Updated SMARTS pattern to better represent the core structure of rotenoids.
    # Taking into account the cis-fused aromatic rings and chromene features.
    # The pattern is a guess based on common rotenoid substructures and may need adjustments
    rotenoid_pattern = Chem.MolFromSmarts("c1ccc2c(c1)C(=O)c3cc4ccc(OC)c4c5c(cc3O2)c(c5)O")

    if rotenoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Match against the SMARTS pattern
    if not mol.HasSubstructMatch(rotenoid_pattern):
        return False, "No tetracyclic rotenoid core structure found"
    
    return True, "Contains the tetracyclic rotenoid core structure"