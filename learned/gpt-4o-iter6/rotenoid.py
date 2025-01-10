"""
Classifies: CHEBI:71543 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Define rotenoid skeleton SMARTS pattern
    # This is a simplification due to complexity; real analysis may need detailed work
    rotenoid_pattern = Chem.MolFromSmarts("C1COc2cc3CCCCOc3cc2C1") # Simplified pattern

    # Check if the structure matches the rotenoid skeleton
    if not mol.HasSubstructMatch(rotenoid_pattern):
        return False, "No tetrahydrochromeno[3,4-b]chromene skeleton found"

    # Further check for potential derived substitutions
    # For simplification purposes, detailed substituent analysis is omitted
    # In practical applications, detailed structural analysis is necessary
    # This might include checking for allowable substitutions at specified positions

    return True, "Contains the tetrahydrochromeno[3,4-b]chromene skeleton characteristic of rotenoids"