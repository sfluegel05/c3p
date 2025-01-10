"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    A rotenoid is a member of the class of tetrahydrochromenochromene that 
    consists of a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

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

    # Update the SMARTS pattern to better match the rotenoid features
    # Attempt to emphasize cis-fusion and tetrahydrochromeno[3,4-b]chromene framework
    # Note: Patterns may vary based on how precise inclusion/exclusion criteria become
    rotenoid_pattern = Chem.MolFromSmarts("C1(C(=O)c2cc3occc3c4cc5occc5c24)COC6OC1")

    if rotenoid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Match against the SMARTS pattern
    if not mol.HasSubstructMatch(rotenoid_pattern):
        return False, "No tetrahydrochromeno[3,4-b]chromene skeleton found"
    
    # Further checks: ensuring any features specific to cis-fusion or substitutions
    # could be detailed here, but often these require empirical data

    return True, "Contains the tetrahydrochromeno[3,4-b]chromene skeleton characteristic of rotenoids"

# Consider testing this with known positivenegatives from a curated dataset or known examples