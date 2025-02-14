"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid contains side chains capable of forming hydrogen bonds.

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

    # Generic amino acid backbone pattern: includes amino group, alpha carbon, and carboxyl group
    aa_backbone = Chem.MolFromSmarts("[NX3][C@@H]([*])[CX3](=O)[OX1]")
    if not mol.HasSubstructMatch(aa_backbone):
        return False, "No amino acid backbone found"

    # Polar side chain functional groups
    polar_group_patterns = [
        "[OX2H]",      # Alcohol or hydroxyl group
        "[#7]",        # Nitrogen (amines and amides)
        "[SX2]",       # Sulfhydryl/thiol group
        "[NX2][NX3]",  # Imidazole (as in histidine)
        "C(=O)O",      # Carboxylic acid
        "C(=O)N"       # Amide
    ]

    # Check if side chain contains any polar functional groups
    for pattern in polar_group_patterns:
        polar_group = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(polar_group):
            return True, f"Polar group '{pattern}' found in side chain"

    return False, "No polar groups found in side chain"