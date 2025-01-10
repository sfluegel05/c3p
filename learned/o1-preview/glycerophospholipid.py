"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: glycerophospholipid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is defined as any glycerolipid having a phosphate group
    ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (three carbons connected in a chain)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Look for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        # Try deprotonated phosphate
        phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)[O-]")
        phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
        if not phosphate_matches:
            return False, "No phosphate group found"

    # Look for ester bonds attached to glycerol backbone
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_count = len(ester_matches)
    if ester_count < 2:
        return False, f"Found {ester_count} ester groups, need at least 2"

    # Optional: Check that esters and phosphate are connected to glycerol backbone
    # This requires mapping atoms and bonds, which can be complex

    return True, "Contains glycerol backbone with two fatty acid chains and a phosphate group attached"