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

    # Look for glycerol backbone pattern (C-C-C with exactly 2 OH and 1 attachment)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4;!H0]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
       glycerol_pattern = Chem.MolFromSmarts("[CH2X4;!H0]([OX2H1])[CHX4]([OX2H1])[CH2X4]")
       glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
         glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4]([OX2H1])[CH2X4]([OX2H1])")
         glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"


    # Define the attachment pattern: ester, ether, vinyl ether, or alkyl.
    attachment_pattern = Chem.MolFromSmarts("[CH2X4;!H0]([OX2H1])[CHX4]([OX2H1])[CH2X4]([OX2;!H0]([CX4](=[OX1])~[!#1])|[OX2]([CX4]~[!#1])|[CX3]=[CX3][OX2][CX4]~[!#1]|[CX4]~[!#1])")

    attachment_matches = mol.GetSubstructMatches(attachment_pattern)
    if not attachment_matches:
        attachment_pattern = Chem.MolFromSmarts("[CH2X4;!H0]([OX2H1])[CHX4]([OX2H1])[CH2X4;!H0]([CX4]~[!#1])")
        attachment_matches = mol.GetSubstructMatches(attachment_pattern)
    if not attachment_matches:
        attachment_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4;!H0]([OX2;!H0]([CX4](=[OX1])~[!#1])|[OX2]([CX4]~[!#1])|[CX3]=[CX3][OX2][CX4]~[!#1])")
        attachment_matches = mol.GetSubstructMatches(attachment_pattern)
    if not attachment_matches:
         attachment_pattern = Chem.MolFromSmarts("[CH2X4]([OX2H1])[CHX4]([OX2H1])[CH2X4;!H0]([CX4]~[!#1])")
         attachment_matches = mol.GetSubstructMatches(attachment_pattern)
    
    if len(attachment_matches) != 1:
        return False, f"Found {len(attachment_matches)} attachment points, require exactly 1"

    return True, "Contains a glycerol backbone with a single acyl, alkyl, or alk-1-enyl substituent"