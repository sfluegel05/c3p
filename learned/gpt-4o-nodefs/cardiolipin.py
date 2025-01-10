"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    Cardiolipin is a diphosphatidylglycerol with four acyl chains.

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

    # Check for two phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("O=P([O-])([O-])O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 2"

    # Check for glycerol backbone (C-C-C)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) < 2:
        return False, f"Missing glycerol backbone, found {len(glycerol_matches)}"

    # Check for 4 ester linkages (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 4"

    return True, "Contains structural motifs of cardiolipin (two phosphates, glycerol backbone, four acyl chains)"