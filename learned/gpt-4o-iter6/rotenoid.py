"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Attempt a new, more detailed SMARTS pattern based on ChEMBL data and known features
    # Identify the main rotenoid skeleton, incorporating cis-fusion and specific ring structures.
    # Note: This might not be a complete pattern; additional characterization in real scenarios is necessary
    # This focus attempts to cover more known skeleton patterns found in known rotenoid structures.
    rotenoid_pattern = Chem.MolFromSmarts("O1Cc2c(c1-c1cc(ccc1)c1oc3c(c2)cccc3o1)C(=O)C") # Skeleton pattern example

    # Check if the structure matches the rotenoid core skeleton
    if not mol.HasSubstructMatch(rotenoid_pattern):
        return False, "No tetrahydrochromeno[3,4-b]chromene skeleton found"
    
    # Further analysis could be done if additional molecular features are needed

    return True, "Contains the tetrahydrochromeno[3,4-b]chromene skeleton characteristic of rotenoids"