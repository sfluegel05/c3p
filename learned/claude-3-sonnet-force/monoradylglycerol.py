"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: CHEBI:36777 monoradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol is a glycerol bearing a single acyl, alkyl or alk-1-enyl substituent at an unspecified position.

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

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_match = mol.GetSubstructMatch(glycerol_pattern)
    if not glycerol_match:
        return False, "No glycerol backbone found"
    
    # Look for single acyl, alkyl or alk-1-enyl substituent pattern
    substituent_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])|[OX2][CX4]|[CX3]=[CX3]")
    substituent_matches = mol.GetSubstructMatches(substituent_pattern)
    if len(substituent_matches) != 1:
        return False, f"Found {len(substituent_matches)} substituents, need exactly 1"
    
    # Check that the substituent is attached to the glycerol backbone
    for idx in substituent_matches[0]:
        if idx in glycerol_match:
            return True, "Contains glycerol backbone with a single acyl, alkyl or alk-1-enyl substituent"
    
    return False, "Substituent not attached to glycerol backbone"