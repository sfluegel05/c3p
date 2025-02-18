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

    # Define SMARTS pattern for 1-monoglyceride core structure:
    # - Glycerol backbone (HO-CH2-OH groups)
    # - Exactly one ester group at position 1 (O-C=O)
    # - Two hydroxyl groups at positions 2 and 3
    # [CH2] with ester group -> [CH2]([OX2][C](=O))...
    # Middle CH with hydroxyl -> [CH]([OH])
    # Terminal CH2 with hydroxyl -> [CH2]([OH])
    glyceride_pattern = Chem.MolFromSmarts("[CH2]([OX2][C](=O)*)[CH]([OH])[CH2]([OH])")
    
    # Find matches for the core structure
    matches = mol.GetSubstructMatches(glyceride_pattern)
    
    if not matches:
        # Check alternative orientation where positions might be reversed
        # Some SMILES representations might have different atom order
        alt_pattern = Chem.MolFromSmarts("[CH2]([OH])[CH]([OH])[CH2]([OX2][C](=O)*)")
        alt_matches = mol.GetSubstructMatches(alt_pattern)
        if len(alt_matches) == 0:
            return False, "Glycerol backbone with ester at position 1 not found"
        else:
            return True, "1-monoglyceride with acyl group at position 1 (alternative orientation)"
    
    # Verify exactly one ester is attached to the glycerol backbone
    # (handled by SMARTS pattern matching)
    return True, "1-monoglyceride with acyl group at position 1"