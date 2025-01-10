"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: monoradylglycerol
A lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Look for specific glycerol backbone patterns with substituents
    # Pattern 1: -O-CH2-CH(OH)-CH2-OH (1-substituted)
    # Pattern 2: -O-CH2-CH(OH)-CH2-OH (2-substituted)
    # Pattern 3: -O-CH2-CH(OH)-CH2-OH (3-substituted)
    glycerol_patterns = [
        "[OX2][CH2][CH]([OH1])[CH2][OH1]",  # Generic pattern
        "[OH1][CH2][CH]([OH1])[CH2][OX2]",  # Alternative arrangement
    ]
    
    found_glycerol = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_glycerol = True
            break
            
    if not found_glycerol:
        return False, "No glycerol backbone with correct substitution pattern found"

    # Count free hydroxyls (-OH groups)
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches != 2:
        return False, f"Found {oh_matches} free hydroxyl groups, need exactly 2"

    # Look for valid substituents
    # Acyl group (ester)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    # Alkyl group (ether)
    ether_pattern = Chem.MolFromSmarts("[OX2]([CH2][CH]([OH1])[CH2][OH1])[CX4]")
    # Alk-1-enyl group
    alkenyl_pattern = Chem.MolFromSmarts("[OX2]([CH2][CH]([OH1])[CH2][OH1])[CH2][CH]=[CH]")
    
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    alkenyl_matches = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    # Must have exactly one type of substituent
    if ester_matches > 0:
        if ester_matches == 1 and ether_matches == 0 and alkenyl_matches == 0:
            return True, "Valid monoradylglycerol with one acyl substituent"
    elif ether_matches > 0:
        if ether_matches == 1 and alkenyl_matches == 0:
            return True, "Valid monoradylglycerol with one alkyl substituent"
    elif alkenyl_matches == 1:
        return True, "Valid monoradylglycerol with one alk-1-enyl substituent"
    
    # Check for carboxylic acids (to avoid false positives)
    acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains carboxylic acid group, not a monoradylglycerol"
        
    # If we get here, the substitution pattern is wrong
    return False, "Does not have correct single substituent pattern"