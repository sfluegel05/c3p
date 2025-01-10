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
    A thiol has an -SH group attached to a carbon atom.

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

    # Look for -SH groups connected to carbon (broader pattern)
    thiol_pattern = Chem.MolFromSmarts("[#6;!$(C(=O)S)]-[SX2H1]")
    
    # Find all matches
    matches = mol.GetSubstructMatches(thiol_pattern)
    
    if not matches:
        return False, "No thiol group (-SH) found attached to carbon"
    
    # Additional checks to exclude false positives
    
    # Check for peptide patterns (to exclude cysteine in peptides)
    peptide_pattern = Chem.MolFromSmarts("[NX3H1,NX3H2][CX4H1][CX3](=[OX1])[NX3H1,NX3H2]")
    cysteine_pattern = Chem.MolFromSmarts("[NX3H1,NX3H2][CX4H1](CS)[CX3](=[OX1])[OH,O-,NX3H1,NX3H2]")
    
    if mol.HasSubstructMatch(peptide_pattern) and mol.HasSubstructMatch(cysteine_pattern):
        return False, "Appears to be a cysteine-containing peptide rather than a free thiol"
    
    # Check for disulfides - these should not be counted as thiols
    disulfide_pattern = Chem.MolFromSmarts("[SX2][SX2]")
    if mol.HasSubstructMatch(disulfide_pattern):
        disulfide_matches = mol.GetSubstructMatches(disulfide_pattern)
        thiol_s_atoms = set(match[1] for match in matches)
        disulfide_s_atoms = set(atom for match in disulfide_matches for atom in match)
        if thiol_s_atoms & disulfide_s_atoms:
            return False, "Contains disulfide bonds rather than free thiol"
    
    # Check that sulfur is not part of a thioester
    thioester_pattern = Chem.MolFromSmarts("[S][CX3](=O)[#6]")
    if mol.HasSubstructMatch(thioester_pattern):
        thioester_matches = mol.GetSubstructMatches(thioester_pattern)
        thiol_s_atoms = set(match[1] for match in matches)
        thioester_s_atoms = set(match[0] for match in thioester_matches)
        if thiol_s_atoms & thioester_s_atoms:
            return False, "Contains thioester rather than free thiol"
            
    # Check for thioacetals and similar structures
    thioacetal_pattern = Chem.MolFromSmarts("[#6]-[SX2]-[#6;!$(C[SX2H1])]")
    if mol.HasSubstructMatch(thioacetal_pattern):
        thioacetal_matches = mol.GetSubstructMatches(thioacetal_pattern)
        thiol_s_atoms = set(match[1] for match in matches)
        thioacetal_s_atoms = set(match[1] for match in thioacetal_matches)
        if thiol_s_atoms & thioacetal_s_atoms:
            return False, "Contains thioacetal rather than free thiol"
    
    # Count the number of thiol groups
    num_thiols = len(matches)
    
    # Success message includes count if more than one thiol group
    if num_thiols > 1:
        return True, f"Contains {num_thiols} thiol groups (-SH) attached to carbon atoms"
    return True, "Contains a thiol group (-SH) attached to a carbon atom"