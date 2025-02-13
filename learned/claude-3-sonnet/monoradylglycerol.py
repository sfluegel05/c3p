"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
A lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (C-C-C with hydroxyls)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Count free hydroxyls (-OH groups)
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches != 2:
        return False, f"Found {oh_matches} free hydroxyl groups, need exactly 2"
    
    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ether_pattern = Chem.MolFromSmarts("[OX2]([CX4])[CX4]")  # For alkyl substituents
    alkenyl_pattern = Chem.MolFromSmarts("[OX2][CH2][CH]=[CH]")  # For alk-1-enyl substituents
    
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    alkenyl_matches = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    total_substituents = ester_matches + ether_matches + alkenyl_matches
    
    # Must have exactly one substituent
    if total_substituents == 0:
        return False, "No acyl, alkyl, or alk-1-enyl substituent found"
    elif total_substituents > 1:
        return False, f"Found {total_substituents} substituents, need exactly 1"
    
    # Verify carbon count (at least 5 for glycerol + smallest possible substituent)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 5:
        return False, "Too few carbons for monoradylglycerol"
    
    # Verify oxygen count (at least 3 for glycerol, plus 1-2 for substituent)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3 or o_count > 5:
        return False, f"Invalid oxygen count ({o_count}) for monoradylglycerol"
    
    # Success - determine type of substituent
    if ester_matches == 1:
        subst_type = "acyl"
    elif ether_matches == 1:
        subst_type = "alkyl"
    else:
        subst_type = "alk-1-enyl"
        
    return True, f"Valid monoradylglycerol with one {subst_type} substituent"