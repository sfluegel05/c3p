"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have side chains that can form hydrogen bonds, such as
    hydroxyl, amides, carboxyl, or basic nitrogen groups.

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

        # Identify amino acid backbone: Look for NH2-C-C(=O)-O
        backbone_pattern = Chem.MolFromSmarts("[NX3H2][CX4][CX3](=O)[OX1H,OX2]")
        if not mol.HasSubstructMatch(backbone_pattern):
            return False, "No amino acid backbone found"

        # Identify polar side chains capable of forming hydrogen bonds
        polar_patterns = [
            Chem.MolFromSmarts("[OX2H]"),  # hydroxyl group
            Chem.MolFromSmarts("[CX3](=O)[NX3]"),  # amide group
            Chem.MolFromSmarts("[CX3](=O)[OX2H]"),  # carboxyl group in R group
            Chem.MolFromSmarts("[NX3][CX3](=O)"),  # basic nitrogen, e.g., arginine
            Chem.MolFromSmarts("[nX2]"),  # heterocyclic nitrogen, e.g., in histidine
            Chem.MolFromSmarts("[SX2H]")  # thiol group, e.g., in cysteine
        ]

        for pattern in polar_patterns:
            if mol.HasSubstructMatch(pattern):
                return True, "Contains a polar side chain capable of forming hydrogen bonds"

        return False, "No polar side chain found capable of forming hydrogen bonds"
    
    except Exception as e:
        return None, f"Error in processing SMILES: {str(e)}"