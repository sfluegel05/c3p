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

    # SMARTS pattern for glycerol backbone with ester at position 1
    # [CH2] (connected to ester O-C=O) - [CH] (OH) - [CH2] (OH)
    glyceride_pattern = Chem.MolFromSmarts(
        "[CH2]([OX2][C](=O))-[CH]([OH])-[CH2]([OH])"
    )
    
    # Check for the pattern
    matches = mol.GetSubstructMatches(glyceride_pattern)
    if not matches:
        return False, "Glycerol backbone with ester at position 1 not found"

    # Check that there's exactly one ester group attached to the glycerol
    # Each match corresponds to a glycerol with one ester, so multiple matches imply multiple glycerol backbones (unlikely)
    if len(matches) > 1:
        return False, "Multiple glycerol backbones found"
    
    # Check substituents on glycerol carbons
    # Get the three carbons in the match
    backbone_indices = matches[0][:3]  # Assuming the first three atoms are the glycerol carbons
    positions = {
        0: "1st (ester)",
        1: "2nd (OH)",
        2: "3rd (OH)"
    }
    
    for i, pos in enumerate(backbone_indices):
        atom = mol.GetAtomWithIdx(pos)
        # Get neighbors not part of the backbone
        neighbors = [n for n in atom.GetNeighbors() if n.GetIdx() not in backbone_indices]
        if i == 0:  # Position 1: should have exactly one ester group
            if len(neighbors) != 1:
                return False, f"Position 1 has {len(neighbors)} substituents (expected 1 ester)"
            ester_oxygen = neighbors[0]
            if ester_oxygen.GetAtomicNum() != 8 or not ester_oxygen.GetBonds()[0].GetBondType() == Chem.BondType.SINGLE:
                return False, "Position 1 substituent is not an ester oxygen"
        elif i == 1:  # Position 2: should have one hydroxyl
            if len(neighbors) != 1 or neighbors[0].GetAtomicNum() != 8:
                return False, "Position 2 has invalid substituent (expected OH)"
        elif i == 2:  # Position 3: should have one hydroxyl
            if len(neighbors) != 1 or neighbors[0].GetAtomicNum() != 8:
                return False, "Position 3 has invalid substituent (expected OH)"

    return True, "1-monoglyceride with acyl group at position 1"