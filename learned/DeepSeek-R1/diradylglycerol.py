"""
Classifies: CHEBI:76578 diradylglycerol
"""
"""
Classifies: CHEBI:64435 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol has a glycerol backbone with exactly two substituent groups
    (acyl, alkyl, or alk-1-enyl) attached via oxygen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Adjusted glycerol pattern to allow 3 bonds on central carbon
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX3][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Get all possible glycerol backbones (may have multiple conformations)
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone found"

    # Define substituent patterns (acyl, alkyl, alk-1-enyl)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3]=[OX1]")  # Acyl (ester)
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")        # Alkyl (ether)
    vinyl_ether_pattern = Chem.MolFromSmarts("[OX2]C=C")    # Alk-1-enyl (vinyl ether)

    # Check all possible glycerol backbones
    for backbone in matches:
        substituent_count = 0
        # Check each carbon in the glycerol backbone
        for carbon_idx in backbone:
            atom = mol.GetAtomWithIdx(carbon_idx)
            # Look for oxygen substituents
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    o_idx = neighbor.GetIdx()
                    # Check if oxygen is part of any substituent pattern
                    if any(o_idx in match for match in mol.GetSubstructMatches(ester_pattern)):
                        substituent_count += 1
                    elif any(o_idx in match for match in mol.GetSubstructMatches(ether_pattern)):
                        substituent_count += 1
                    elif any(o_idx in match for match in mol.GetSubstructMatches(vinyl_ether_pattern)):
                        substituent_count += 1

        # Check for exactly two substituents in this backbone configuration
        if substituent_count == 2:
            return True, "Contains glycerol backbone with exactly two substituent groups (acyl/alkyl/alk-1-enyl)"
    
    return False, f"Found {substituent_count} substituent groups, need exactly 2"