"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid usually has a cyclopenta[a]phenanthrene carbon skeleton,
    combined with typical substitutions and stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid SMARTS: 3 hexane rings and 1 pentane ring (start with a core structure)
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C4CCC(C4)C3C2C1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain cyclopenta[a]phenanthrene skeleton"
    
    # Check for at least some stereo centers. Steroids typically have quite a few
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 2:
        return False, "Not enough stereochemistry typically found in steroids"
    
    # Check for common functional groups attached to steroids (e.g., hydroxyl, keto)
    functional_group_patterns = [
        Chem.MolFromSmarts("[OX2H]"),  # hydroxyl groups
        Chem.MolFromSmarts("[CX3](=O)"),  # keto groups
    ]
    for pattern in functional_group_patterns:
        if mol.HasSubstructMatch(pattern):
            break
    else:
        return False, "Missing typical functional groups found in steroids (e.g. hydroxyl, keto)"

    return True, "Contains features typical of a steroid: cyclopenta[a]phenanthrene skeleton with appropriate stereochemistry and functional groups"