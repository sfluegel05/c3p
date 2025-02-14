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
    
    # Define glycerol backbone pattern with phosphate group ester-linked to terminal carbon
    glycerol_phosphate_pattern = Chem.MolFromSmarts("""
    [C@@H]([O]*)([C@@H]([O]*)([C@@H]([O]*)([O][P](=O)([O-])[O*])))
    """)
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Glycerol backbone with phosphate group not found"
    
    # Check for ester or ether linkages at positions 1 and 2
    # Position 1 linkage (can be ester or ether)
    pos1_ester = Chem.MolFromSmarts('[O][CH2][C](=O)[C]')
    pos1_ether = Chem.MolFromSmarts('[O][CH2][CH][C]')
    pos1_match = mol.HasSubstructMatch(pos1_ester) or mol.HasSubstructMatch(pos1_ether)
    if not pos1_match:
        return False, "No ester or ether linkage found at position 1 of glycerol backbone"

    # Position 2 linkage (can be ester or ether)
    pos2_ester = Chem.MolFromSmarts('[O][CH][C](=O)[C]')
    pos2_ether = Chem.MolFromSmarts('[O][CH][CH][C]')
    pos2_match = mol.HasSubstructMatch(pos2_ester) or mol.HasSubstructMatch(pos2_ether)
    if not pos2_match:
        return False, "No ester or ether linkage found at position 2 of glycerol backbone"

    # Check that phosphate group is ester-linked to terminal carbon
    phosphate_pattern = Chem.MolFromSmarts('[O][CH2][O][P](=O)([O-])[O*]')
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group ester-linked to terminal carbon not found"

    return True, "Contains glycerol backbone with phosphate group ester-linked to terminal carbon"