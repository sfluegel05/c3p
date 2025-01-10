"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is defined as a monoglyceride with an acyl group attached to the primary hydroxyl group of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1-monoglyceride, False otherwise
        str: Explanation of the result
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern: primary alcohols adjacent to secondary alcohol
    glycerol_backbone = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_backbone):
        return False, "No glycerol backbone found"

    # Check for ester group at the primary hydroxyl position
    ester_primary_position = Chem.MolFromSmarts("COC(=O)[CX4H]")
    ester_matches = mol.GetSubstructMatches(ester_primary_position)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1 attached to primary alcohol"

    # Check for stereochemistry if specified as sn-glycerol
    chiral_center_check = any(atom.GetChiralTag() for atom in mol.GetAtoms())
    sn_pattern = Chem.MolFromSmarts("[C@H](O)C(O)")
    sn_match = mol.HasSubstructMatch(sn_pattern)
    
    if chiral_center_check and not sn_match:
        return False, "Found stereochemistry but not the expected sn configuration"

    return True, "Contains glycerol backbone with an acyl group esterified at the primary position, consistent with a 1-monoglyceride"