"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have a carboxyl and amine functional group, with side chains that can form hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alpha-amino acid core structure: N-C-C(=O)-O
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=O)[OX1]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No alpha-amino acid core structure found"

    # Check for polar side chains: OH, CONH2, SH, NH2
    polar_side_chain_patterns = [
        Chem.MolFromSmarts("[OX2H]"),   # OH group
        Chem.MolFromSmarts("[CX3](=O)[NH2]"),  # Amide group
        Chem.MolFromSmarts("[SX2H]"),   # SH group
        Chem.MolFromSmarts("[NX3H2]")   # NH2 group
    ]
    
    for pattern in polar_side_chain_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains polar side chain capable of hydrogen bonding"

    return False, "No polar side chain found capable of hydrogen bonding"