"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have side chains that can form hydrogen bonds, such as hydroxyl, amides, carboxyl, or basic nitrogen groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
        
        # Identify amino acid backbone (more flexible to allow variations)
        backbone_pattern = Chem.MolFromSmarts("[NX3][CH1]([C;R0])[C](=O)O")  # N-C-(C)-C(=O)-O pattern without wrapping, R0 means non-ring
        if not mol.HasSubstructMatch(backbone_pattern):
            return False, "No amino acid backbone found"

        # Identify polar side chain groups
        polar_patterns = [
            Chem.MolFromSmarts("[$([NX3][CX3](=O)[OX1]),$([CX3](=O)[OH]),$([OH]),$([NX3H]),$([nX2]) ,$([N+]),$([SX2H])]")  # Extended patterns
        ]

        # Check for any polar pattern match
        for pattern in polar_patterns:
            if mol.HasSubstructMatch(pattern):
                return True, "Contains a polar side chain capable of forming hydrogen bonds"

        return False, "No polar side chain found capable of forming hydrogen bonds"
    
    except Exception as e:
        return None, f"Error in processing SMILES: {str(e)}"