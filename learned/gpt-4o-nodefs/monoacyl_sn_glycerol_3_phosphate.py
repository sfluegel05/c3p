"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a glycerol backbone with the correct stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO)[C@@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No stereo-specific glycerol backbone found"

    # Correct acyl chain detection with ester bond C(=O)O
    acyl_pattern = Chem.MolFromSmarts("[C,C@H](=O)OC")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) != 1:
        return False, f"Expected 1 acyl chain, found {len(acyl_matches)}"

    # Check phosphate group connection at the third position
    phosphate_pattern = Chem.MolFromSmarts("[O]P(=O)(O)OC")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found attached to glycerol at the correct position"
    
    return True, "Contains stereo-specific glycerol, single acyl chain, and phosphate group at the third position"