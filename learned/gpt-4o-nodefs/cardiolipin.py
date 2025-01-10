"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    Cardiolipin consists of a diphosphatidylglycerol structure with four acyl chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for two phosphate groups, consider various representations
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(O)(=O)O"), # Neutral phosphate
        Chem.MolFromSmarts("O=P([O-])(O)O"), # Anionic/diester phosphate
        Chem.MolFromSmarts("O=P(O)(O)O")  # General phosphate representation
    ]
    phosphate_count = sum(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns)
    if phosphate_count != 2:
        return False, f"Found {phosphate_count} phosphate groups, need exactly 2"

    # Check for the general glycerol backbone pattern, considering stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](CO)O")  # Include a central stereochemistry
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) < 1:
        return False, f"Missing or incorrect glycerol backbone pattern, found {len(glycerol_matches)}"

    # Check typically for 4 ester linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 4"

    return True, "Contains structural motifs of cardiolipin (two phosphates, glycerol backbone, four acyl chains)"