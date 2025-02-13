"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride
Definition: A monoglyceride in which the acyl substituent is located at position 1
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has:
    - A glycerol backbone (3 carbons)
    - An ester group at position 1
    - Hydroxyl groups at positions 2 and 3
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_1_monoglyceride, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific pattern for 1-monoglyceride:
    # [CH2] with ester at position 1
    # connected to [CH] with OH
    # connected to [CH2] with OH
    # The pattern ensures correct connectivity and substitution
    mono_1_pattern = Chem.MolFromSmarts("""
        [CH2X4]-[OX2]-[CX3](=[OX1])-[*] # Position 1 with ester
        .[CH2X4]-[OX2H1] # Position 3 with OH
        .[CHX4]-[OX2H1] # Position 2 with OH
        """)
    
    if not mol.HasSubstructMatch(mono_1_pattern):
        return False, "Does not match 1-monoglyceride pattern"

    # Verify it's a single glycerol backbone by checking connectivity
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) != 1:
        return False, "Must have exactly one glycerol backbone"

    # Count ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Verify the ester is connected to a carbon chain (not just a small group)
    acyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4,CX3]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "Ester must be connected to a carbon chain"

    # Additional check to ensure the glycerol backbone is properly connected
    backbone_pattern = Chem.MolFromSmarts("""
        [CH2X4]-[OX2]-[CX3](=[OX1]) # Position 1
        -[CHX4]-[OX2H1] # Position 2
        -[CH2X4]-[OX2H1] # Position 3
        """)
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Glycerol backbone not properly connected"

    return True, "Contains glycerol backbone with one ester bond at position 1 and two free hydroxyl groups"