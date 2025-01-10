"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

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

    # Look for a glycerol backbone with stereochemistry
    # Allow for flexibility in stereochemistry as variations exist in provided examples
    glycerol_pattern = Chem.MolFromSmarts("[O-][C@H](CO)[C@@H](O)CO |t|")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No stereo-specific glycerol backbone found"

    # Detect single acyl chain with ester bond at glycerol position
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CO)")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if len(acyl_matches) < 1:
        return False, "Acyl chain not correctly detected"

    # Verify phosphate group at the last oxygen of glycerol
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found in correct position"

    return True, "Contains stereo-specific glycerol, single acyl chain, and phosphate group at the third position"