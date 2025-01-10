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
    
    # Improved glycerol-like backbone pattern including stereochemistry
    glycerol_like_pattern = Chem.MolFromSmarts("[C@@H]([OH1])-[C@H]([OH1])-O")
    if not mol.HasSubstructMatch(glycerol_like_pattern):
        return False, "No suitable glycerol-like backbone found"

    # Look for ether linkage pattern - R-O-R
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkage found"

    # Checking for phosphate group is optional, pay attention to its presence
    phospho_pattern = Chem.MolFromSmarts("P(=O)([O-])([O-])")
    if mol.HasSubstructMatch(phospho_pattern):
        phosphate_presence = " and phosphate group identified"
    else:
        phosphate_presence = ""

    # If there are more ether than ester linkages, it's likely an ether lipid
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ether_matches) > len(ester_matches):
        return True, f"Ether lipid identified by ether linkages on glycerol-like backbone{phosphate_presence}"
    else:
        return False, "Ether-to-ester linkage ratio not sufficient for ether lipid classification"