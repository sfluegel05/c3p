"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    Phosphatidylinositol phosphates are characterized by a glycerol backbone with two fatty acid chains,
    a phosphate group, and an inositol ring that is phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is phosphatidylinositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General check for glycerol backbone with ester linkages
    glycerol_pattern = Chem.MolFromSmarts("COC([O])([C@H])CO")  # Broad pattern to include variations
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No complete glycerol backbone with proper ester linkages found"
        
    # Look for fatty acid ester groups (More flexibility with stereo and bonds)
    ester_pattern = Chem.MolFromSmarts("OC(=O)C")  # General ester pattern
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Inositol and phosphate group structure recognition
    # Updated pattern to recognize varied stereo and phosphated positions
    inositol_phosphate_pattern = Chem.MolFromSmarts("C1([O-,O])C([O-,O])C([O-,O])C([O-,O])C([O-,O])C1O[PX4](=[OX1])([OX2])[OX1-]")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "Inositol ring phosphate structure incomplete or missing"

    # Check for phosphate groups: more than one can be present
    phosphate_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("P(=O)([O-,O])[O-,O]")))
    if phosphate_count < 1:
        return False, "Insufficient number of phosphate groups (need at least one)"

    return True, "Contains features characteristic of phosphatidylinositol phosphate"