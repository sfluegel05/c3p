"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino acid based on its SMILES string.
    An N-hydroxy-alpha-amino acid is an amino acid in which at least one hydrogen 
    attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Combined SMARTS pattern for N-hydroxy-alpha-amino acids, including both -N(OH)- and -N(OH)2 forms
    # The key here is to ensure that the -N(OH) or -N(OH)2 is directly connected to the alpha carbon of the amino acid.
    # -N(OH)- and -N(OH)2 directly connected to the alpha carbon and a carbonyl group (C=O)
    n_hydroxy_alpha_amino_acid_pattern1 = Chem.MolFromSmarts("[NX3]([OX2])([H])-[CX4]([H])(-[CX3](=[OX1])-[OX2])")
    n_hydroxy_alpha_amino_acid_pattern2 = Chem.MolFromSmarts("[NX3]([OX2])([OX2])-[CX4]([H])(-[CX3](=[OX1])-[OX2])")
    # Tautomers where the N is protonated 
    n_hydroxy_alpha_amino_acid_pattern3 = Chem.MolFromSmarts("[NX3+]([OX2][H])([H])-[CX4]([H])(-[CX3](=[OX1])-[OX2])")
    n_hydroxy_alpha_amino_acid_pattern4 = Chem.MolFromSmarts("[NX3+]([OX2][H])([OX2][H])-[CX4]([H])(-[CX3](=[OX1])-[OX2])")
    if mol.HasSubstructMatch(n_hydroxy_alpha_amino_acid_pattern1) or \
       mol.HasSubstructMatch(n_hydroxy_alpha_amino_acid_pattern2) or \
       mol.HasSubstructMatch(n_hydroxy_alpha_amino_acid_pattern3) or \
       mol.HasSubstructMatch(n_hydroxy_alpha_amino_acid_pattern4):
        return True, "Contains an alpha amino acid structure with at least one hydroxy group on the amino nitrogen"
    
    return False, "No N-hydroxy group attached to the alpha amino group found"