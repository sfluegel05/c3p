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

    # Refined steroid pattern: a more generalized cyclopenta[a]phenanthrene skeleton
    steroid_pattern = Chem.MolFromSmarts("C1C2C3C4C5=C(C3C2)C2CC(C1)CCC2CCC4")  # more flexible for steroids
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain cyclopenta[a]phenanthrene skeleton"
    
    # Check for some stereo centers, but not too strictly
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 2:
        return False, "Not enough stereochemistry typically found in steroids"

    # Broaden check for functional groups typical in steroids
    functional_group_patterns = [
        Chem.MolFromSmarts("[OX2H]"),  # hydroxyl groups
        Chem.MolFromSmarts("[CX3](=O)"),  # keto groups
        Chem.MolFromSmarts("[C;R0]=O"),  # potential ester bond
    ]
    for pattern in functional_group_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains features typical of a steroid: cyclopenta[a]phenanthrene skeleton with appropriate stereochemistry and functional groups"

    return False, "Missing typical functional groups found in steroids (e.g., hydroxyl, keto, ester)"