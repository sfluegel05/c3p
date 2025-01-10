"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid typically has a carbon skeleton derived from a monoterpene
    with around 10 carbons, potentially rearranged or modified by the removal of
    one or more skeletal atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for monoterpenoid carbon count range
    if c_count < 8 or c_count > 13:
        return False, f"Carbon count {c_count} not typical for monoterpenoid"

    # Search for common skeletal patterns (e.g., isoprene units)
    isoprene_pattern = Chem.MolFromSmarts("C(C)(C)C=C")  # Basic isoprene unit
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene unit found, not typical for a monoterpenoid"

    # Check for functional groups typically found in monoterpenoids
    # Alcohols
    alcohol_pattern = Chem.MolFromSmarts("[OX2H]")
    if mol.HasSubstructMatch(alcohol_pattern):
        return True, "Contains isoprene unit and alcohol functional group, consistent with monoterpenoid"

    # Ketones
    ketone_pattern = Chem.MolFromSmarts("C(=O)C")
    if mol.HasSubstructMatch(ketone_pattern):
        return True, "Contains isoprene unit and ketone functional group, consistent with monoterpenoid"

    # Esters
    ester_pattern = Chem.MolFromSmarts("COC(=O)")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Contains isoprene unit and ester functional group, consistent with monoterpenoid"

    return False, "Unable or insufficient evidence to classify as monoterpenoid"