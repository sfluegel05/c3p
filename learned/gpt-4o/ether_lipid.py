"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has one or more carbon atoms on the glycerol or similar backbone linked
    to an alkyl chain via an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Revised glycerol-like backbone pattern without stereochemistry
    glycerol_like_pattern = Chem.MolFromSmarts("[CX4][CX4][OX2H1]")
    if not mol.HasSubstructMatch(glycerol_like_pattern):
        return False, "No suitable glycerol-like backbone found"
    
    # Look for ether linkage pattern - more generalized
    ether_pattern = Chem.MolFromSmarts("[CX4]O[CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkage found"

    # Simplified phosphate group detection (not a requirement for all ether lipids)
    phospho_pattern = Chem.MolFromSmarts("[PX4](=O)(O)(O)")
    if mol.HasSubstructMatch(phospho_pattern):
        phosphate_presence = " and phosphate group identified"
    else:
        phosphate_presence = ""

    # Count ester groups to balance classification more broadly
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if len(ether_matches) >= 1 and len(ether_matches) > len(ester_matches):
        return True, f"Ether lipid identified by ether linkages{phosphate_presence}"
    else:
        return False, "Ether linkage presence insufficient for ether lipid classification"