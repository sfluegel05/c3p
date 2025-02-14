"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an amino acid where the amine group is acetylated
    (CH3-C(=O)-N).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for N-acetyl amino acid
    acetyl_amino_acid_pattern = Chem.MolFromSmarts("[NX3](C(=O)C)[CX4][CX3](=[OX1])[OX2H0]")

    if mol.HasSubstructMatch(acetyl_amino_acid_pattern):
        return True, "Molecule matches the N-acetyl-amino acid pattern."
    else:
        return False, "Molecule does not match the N-acetyl-amino acid pattern."