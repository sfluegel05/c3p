"""
Classifies: CHEBI:64611 ether lipid
"""
from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid has one or more carbon atoms on a glycerol or similar backbone linked
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
    
    # Revised glycerol-like backbone pattern, with more flexibility and stereo options
    glycerol_like_pattern = Chem.MolFromSmarts("[O,C]C([O,C])C([O,C])")
    if not mol.HasSubstructMatch(glycerol_like_pattern):
        return False, "No suitable glycerol-like backbone found"

    # Ether linkage detection, considering flexibility in connectivity
    ether_pattern = Chem.MolFromSmarts("[CX4]O[CX4,CX3]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkage found"

    # Simplified phosphate group detection (informative but not required)
    phospho_pattern = Chem.MolFromSmarts("[PX4](=O)(O)(O)(O)")
    phosphate_presence = mol.HasSubstructMatch(phospho_pattern)

    # Determine the ratio of ether versus ester linkage matches
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    if len(ether_matches) >= 1 and len(ether_matches) > len(ester_matches) / 2:
        reason = "Ether lipid identified by ether linkages"
        if phosphate_presence:
            reason += " and phosphate group identified"
        return True, reason
    else:
        return False, "Ether linkage presence insufficient for ether lipid classification"