"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol with a single acyl, alkyl, or alk-1-enyl substituent.

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

    # Look for glycerol backbone pattern (C-C-C with two OH groups)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4]")
    if glycerol_pattern is None:
      return False, "Invalid glycerol SMARTS pattern"
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Define attachment patterns
    ester_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4][OX2][CX3](=[OX1])~[!#1]")
    ether_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4][OX2][CX4]~[!#1]")
    alkyl_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4][CX4]~[!#1]")
    vinyl_ether_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4][OX2][CX3]=[CX3]~[!#1]")


    if ester_pattern is None or ether_pattern is None or alkyl_pattern is None or vinyl_ether_pattern is None:
         return False, "Invalid attachment SMARTS pattern"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    vinyl_ether_matches = mol.GetSubstructMatches(vinyl_ether_pattern)

    total_matches = len(ester_matches) + len(ether_matches) + len(alkyl_matches) + len(vinyl_ether_matches)

    if total_matches != 1:
        return False, f"Found {total_matches} attachment points, require exactly 1"

    return True, "Contains a glycerol backbone with a single acyl, alkyl, or alk-1-enyl substituent"