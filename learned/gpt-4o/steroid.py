"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid usually has a cyclopenta[a]phenanthrene carbon skeleton.

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

    # Improved cyclopenta[a]phenanthrene skeleton: 4 rings (3x6-membered, 1x5-membered)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3CCC4CCCC5=C4C3=C2C1C5")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain cyclopenta[a]phenanthrene skeleton"
    
    # Check for at least some stereochemistry, but not overly strict
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=False)
    if len(chiral_centers) < 2:
        return False, "Not enough stereochemistry typically found in steroids"

    return True, "Contains features typical of a steroid: cyclopenta[a]phenanthrene skeleton"