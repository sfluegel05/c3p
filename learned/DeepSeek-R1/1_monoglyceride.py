"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: 1-monoglyceride (CHEBI:17408)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has a glycerol backbone with a single acyl group at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one ester group (-O-C=O-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check glycerol backbone with ester at position 1 (O-C=O attached to first carbon)
    # SMARTS: [CH2](O-C=O) connected to CH(OH) and CH2(OH)
    glycerol_smarts = Chem.MolFromSmarts("[CH2]([OX2][C](=O)*)[CH]([OH])[CH2]([OH])")
    if not mol.HasSubstructMatch(glycerol_smarts):
        return False, "Glycerol backbone with ester at position 1 not found"

    return True, "1-monoglyceride with acyl group at position 1"