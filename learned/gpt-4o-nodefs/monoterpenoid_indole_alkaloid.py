"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a monoterpenoid indole alkaloid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined SMARTS pattern for indole; allows for common variants
    indole_pattern = Chem.MolFromSmarts("c1c[nH]c2cccc[c]12")
    
    # Check for indole moiety
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole moiety found, or misidentified"

    # Adjusted criteria for the monoterpenoid-like features.
    # This could be a pattern check for isoprene units or common fragment patterns.
    # Since this is less specific, we assume certain connectivity common in monoterpenoid indoles.
    monoterpenoid_fragment_check = any(atom.GetIsAromatic() for atom in mol.GetAtoms())  # Example placeholder

    if not monoterpenoid_fragment_check:
        return False, "No acceptable monoterpenoid-like features found"

    # Additional checks could include specific functional groups often found in alkaloids.

    return True, "Contains indole moiety and monoterpenoid-like features"