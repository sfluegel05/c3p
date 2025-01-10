"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid must have a glycerol backbone with two fatty acid chains
    attached via ester linkages and a phosphate group ester-linked to a terminal carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification or failure
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (O-C-C-C-O)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)COP")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with phosphate found"

    # Look for two ester or ether groups (-O-C(=O)- or -O-C-)
    ester_or_ether_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[CX4]")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    ester_matches = mol.GetSubstructMatches(ester_or_ether_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Ensure there are at least two ester or ether groups
    if len(ester_matches) + len(ether_matches) < 2:
        return False, f"Found {len(ester_matches) + len(ether_matches)} ester or ether groups, need at least 2"

    # Check for phosphate group at terminal position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No terminal phosphate group found"

    return True, "Contains glycerol backbone with two fatty acid chains and terminal phosphate group"