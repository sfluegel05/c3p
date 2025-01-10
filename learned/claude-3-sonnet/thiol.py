"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol has an -SH group attached to a carbon atom of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for -SH groups connected to carbon
    # [#6] is carbon
    # [SX2H1] is -SH group (sulfur with 2 connections, one being H)
    # !$(C(=O)S) excludes thioacids
    thiol_pattern = Chem.MolFromSmarts("[#6;!$(C(=O)S)][SX2H1]")
    
    # Find all matches
    matches = mol.GetSubstructMatches(thiol_pattern)
    
    if not matches:
        return False, "No thiol group (-SH) found attached to carbon"
    
    # Get atom indices of all sulfur atoms that are part of thiol groups
    thiol_s_atoms = set(match[1] for match in matches)
    
    # Check for disulfides - exclude if the thiol sulfur is part of a disulfide
    disulfide_pattern = Chem.MolFromSmarts("[SX2][SX2]")
    if mol.HasSubstructMatch(disulfide_pattern):
        disulfide_matches = mol.GetSubstructMatches(disulfide_pattern)
        disulfide_s_atoms = set(atom for match in disulfide_matches for atom in match)
        if thiol_s_atoms & disulfide_s_atoms:
            return False, "Contains disulfide bonds rather than free thiol"
    
    # Check for thioesters
    thioester_pattern = Chem.MolFromSmarts("[#6][SX2]C(=O)[#6]")
    if mol.HasSubstructMatch(thioester_pattern):
        thioester_matches = mol.GetSubstructMatches(thioester_pattern)
        thioester_s_atoms = set(match[1] for match in thioester_matches)
        if thiol_s_atoms & thioester_s_atoms:
            return False, "Contains thioester rather than free thiol"
    
    # Check for sulfides/thioethers where the sulfur is not a thiol
    thioether_pattern = Chem.MolFromSmarts("[#6][SX2][#6]")
    if mol.HasSubstructMatch(thioether_pattern):
        thioether_matches = mol.GetSubstructMatches(thioether_pattern)
        thioether_s_atoms = set(match[1] for match in thioether_matches)
        # Only exclude if all sulfur atoms are part of thioethers
        if thiol_s_atoms.issubset(thioether_s_atoms):
            return False, "Contains only sulfides/thioethers, no free thiol"
    
    # Count valid thiol groups
    num_thiols = len(matches)
    
    # Success message includes count if more than one thiol group
    if num_thiols > 1:
        return True, f"Contains {num_thiols} thiol groups (-SH) attached to carbon atoms"
    return True, "Contains a thiol group (-SH) attached to a carbon atom"