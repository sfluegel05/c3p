"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have side chains that can form hydrogen bonds, such as hydroxyl, amides, carboxyl, or basic nitrogen groups,
    without forming large peptide or polypeptide structures.

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

        # Define amino acid backbone with limited connectivity to avoid peptides
        backbone_pattern = Chem.MolFromSmarts("N[C@@H]C(=O)O|N[C@H]C(=O)O")  # Address different stereochemistries
        if not mol.HasSubstructMatch(backbone_pattern):
            return False, "No amino acid backbone found"

        # Define polar side chain groups relevant to amino acids
        polar_patterns = [
            Chem.MolFromSmarts("[OX2H]"),    # Hydroxyl (as in serine, tyrosine, threonine)
            Chem.MolFromSmarts("C(=O)[NX3]"),# Amide (as in glutamine, asparagine)
            Chem.MolFromSmarts("C(=O)[OH]"), # Carboxyl on side chain (as in glutamic/aspartic acid)
            Chem.MolFromSmarts("N[C@@H](Cc1ncnc1)[C](=O)"), # His unique polar ring hydrogen bonding
            Chem.MolFromSmarts("C[*:1][NX3H2]")  # Primary amine in side chain (as in lysine)
        ]

        # Ensure polar group is part of the side chain
        for pattern in polar_patterns:
            if mol.HasSubstructMatch(pattern):
                return True, "Contains a polar side chain capable of forming hydrogen bonds"

        return False, "No relevant polar side chain found capable of forming hydrogen bonds"
    
    except Exception as e:
        return None, f"Error in processing SMILES: {str(e)}"