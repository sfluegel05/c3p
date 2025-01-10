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

        # Identify amino acid backbone: Look for NH2-C-COO motif
        backbone_pattern = Chem.MolFromSmarts("[NX3H2,NH,NH3+][CH1](C(=O)O)C")  # N-C-C with a terminal carboxyl group
        if not mol.HasSubstructMatch(backbone_pattern):
            return False, "No amino acid backbone found"

        # Identify polar side chains capable of forming hydrogen bonds
        polar_patterns = [
            Chem.MolFromSmarts("[CX3](=O)[OX2H1]"),  # Carboxylic acids
            Chem.MolFromSmarts("[CX3](=O)[NH2]"),   # Amides
            Chem.MolFromSmarts("[OH1]"),            # Hydroxyl groups
            Chem.MolFromSmarts("[nX2]"),            # Heterocyclic nitrogens, e.g., in histidine
            Chem.MolFromSmarts("[N+H](C)(C)"),      # Charged amine for arginine
            Chem.MolFromSmarts("[SX2H1]"),          # Thiol group, e.g., in cysteine
        ]

        for pattern in polar_patterns:
            if mol.HasSubstructMatch(pattern):
                return True, "Contains a polar side chain capable of forming hydrogen bonds"

        return False, "No polar side chain found capable of forming hydrogen bonds"
    
    except Exception as e:
        return None, f"Error in processing SMILES: {str(e)}"