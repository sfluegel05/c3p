"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: glycerophospholipid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerolipid having a phosphate group ester-linked to a terminal carbon of the glycerol backbone.

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
    
    # Check for exactly one phosphate group
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)(O)')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 1"

    # Define glycerol backbone with phosphate group on a terminal carbon
    glycerol_phosphate_pattern = Chem.MolFromSmarts('[C;H2][C;H][C;H2][O][P](=O)(O)O')
    backbone_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if not backbone_matches:
        # Try the glycerol backbone with phosphate group on the other terminal carbon
        glycerol_phosphate_pattern = Chem.MolFromSmarts('[C;H2][C;H][C;H2][O][C;H2][O][P](=O)(O)O')
        backbone_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
        if not backbone_matches:
            return False, "Glycerol backbone with phosphate group not found"
    
    # Define patterns for ester and ether linkages at positions 1 and 2
    ester_pattern = Chem.MolFromSmarts('[C;H2][O][C](=O)[C]')
    ether_pattern = Chem.MolFromSmarts('[C;H2][O][C;H]')

    # Check for ester or ether linkage at position 1
    position1_matches = mol.GetSubstructMatches(ester_pattern) + mol.GetSubstructMatches(ether_pattern)
    if not position1_matches:
        return False, "No ester or ether linkage found at position 1 of glycerol backbone"

    # Check for ester or ether linkage at position 2
    ester_pattern_pos2 = Chem.MolFromSmarts('[C;H][O][C](=O)[C]')
    ether_pattern_pos2 = Chem.MolFromSmarts('[C;H][O][C;H]')
    position2_matches = mol.GetSubstructMatches(ester_pattern_pos2) + mol.GetSubstructMatches(ether_pattern_pos2)
    if not position2_matches:
        return False, "No ester or ether linkage found at position 2 of glycerol backbone"

    # Ensure that molecule is not a cardiolipin (exclude molecules with more than one phosphate group)
    phosphate_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('P(=O)(O)(O)')))
    if phosphate_count != 1:
        return False, "Molecule has more than one phosphate group, likely not a glycerophospholipid"

    return True, "Contains glycerol backbone with ester or ether-linked fatty acids and a phosphate group at terminal carbon"