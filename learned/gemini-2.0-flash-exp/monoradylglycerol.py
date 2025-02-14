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
    
    # Define glycerol backbone patterns allowing attachment to any of the three carbons
    glycerol_pattern1 = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4][OX2]")
    glycerol_pattern2 = Chem.MolFromSmarts("[CH2X4]([OX2])[CHX4]([OX2H1])[CH2X4]([OX2H1])")
    glycerol_pattern3 = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2])[CH2X4]([OX2H1])")

    if glycerol_pattern1 is None or glycerol_pattern2 is None or glycerol_pattern3 is None:
      return False, "Invalid glycerol SMARTS pattern"
    
    glycerol_matches1 = mol.GetSubstructMatches(glycerol_pattern1)
    glycerol_matches2 = mol.GetSubstructMatches(glycerol_pattern2)
    glycerol_matches3 = mol.GetSubstructMatches(glycerol_pattern3)
    
    if not glycerol_matches1 and not glycerol_matches2 and not glycerol_matches3:
        return False, "No glycerol backbone found"

    # Define attachment patterns
    attachment_pattern1 = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4][OX2][CX3,CX4]~[!#1]")
    attachment_pattern2 = Chem.MolFromSmarts("[CH2X4]([OX2])[CHX4]([OX2H1])[CH2X4]([OX2H1])[CX3,CX4]~[!#1]")
    attachment_pattern3 = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2])[CH2X4]([OX2H1])[CX3,CX4]~[!#1]")

    if attachment_pattern1 is None or attachment_pattern2 is None or attachment_pattern3 is None:
         return False, "Invalid attachment SMARTS pattern"
    
    attachment_matches1 = mol.GetSubstructMatches(attachment_pattern1)
    attachment_matches2 = mol.GetSubstructMatches(attachment_pattern2)
    attachment_matches3 = mol.GetSubstructMatches(attachment_pattern3)
    
    total_matches = len(attachment_matches1) + len(attachment_matches2) + len(attachment_matches3)
    
    if total_matches != 1:
        return False, f"Found {total_matches} attachment points, require exactly 1"

    return True, "Contains a glycerol backbone with a single acyl, alkyl, or alk-1-enyl substituent"