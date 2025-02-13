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

    # Count total oxygens - monoradylglycerols should have 3-4 oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if not (3 <= o_count <= 4):
        return False, f"Invalid number of oxygen atoms ({o_count}), expected 3-4"

    # Look for glycerol backbone patterns
    glycerol_patterns = [
        # 1-substituted
        "[OX2][CH2][CH]([OH1])[CH2][OH1]",
        # 2-substituted
        "[OH1][CH2][CH]([OX2])[CH2][OH1]",
        # 3-substituted
        "[OH1][CH2][CH]([OH1])[CH2][OX2]"
    ]
    
    found_glycerol = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_glycerol = True
            break
            
    if not found_glycerol:
        return False, "No glycerol backbone found"

    # Count free hydroxyls
    oh_pattern = Chem.MolFromSmarts("[OX2H1]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches != 2:
        return False, f"Found {oh_matches} free hydroxyl groups, need exactly 2"

    # Check for valid substituents
    # Acyl group (ester) - more specific pattern
    ester_pattern = Chem.MolFromSmarts("[OX2]([CH2,CH1][CH]([OH1])[CH2][OH1])[CX3](=[OX1])[#6]")
    # Alkyl group (ether) - more specific pattern
    ether_pattern = Chem.MolFromSmarts("[OX2]([CH2,CH1][CH]([OH1])[CH2][OH1])[CX4;!$(C=O)]")
    # Alk-1-enyl group
    alkenyl_pattern = Chem.MolFromSmarts("[OX2]([CH2,CH1][CH]([OH1])[CH2][OH1])[CH2][CH]=[CH]")
    
    ester_matches = len(mol.GetSubstructMatches(ester_pattern))
    ether_matches = len(mol.GetSubstructMatches(ether_pattern))
    alkenyl_matches = len(mol.GetSubstructMatches(alkenyl_pattern))
    
    total_substituents = ester_matches + ether_matches + alkenyl_matches
    
    if total_substituents != 1:
        return False, f"Found {total_substituents} substituents, need exactly 1"

    # Classify based on substituent type
    if ester_matches == 1:
        # Additional check for ester group
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH2,CH1][CH]([OH1])[CH2][OH1]")):
            return True, "Valid monoradylglycerol with one acyl substituent"
    elif ether_matches == 1:
        # Additional check for proper ether linkage
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[OX2]([CH2,CH1][CH]([OH1])[CH2][OH1])[CX4]")):
            return True, "Valid monoradylglycerol with one alkyl substituent"
    elif alkenyl_matches == 1:
        return True, "Valid monoradylglycerol with one alk-1-enyl substituent"

    return False, "Does not have correct single substituent pattern"