"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    This requires an amino acid in which at least one hydrogen attached to the amino group 
    is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # SMARTS pattern improved for alpha amino acid backbone:
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3;!$(N-C=O)]-[C;!H0]-[CX3](=[OX1])-[OX2H1,OX1-]")
    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No alpha amino acid backbone found"
    
    # Search for at least one hydroxy substitution on the N
    n_hydroxy_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]O")
    n_di_hydroxy_pattern = Chem.MolFromSmarts("[NX3](O)O")
    
    n_hydroxy_patterns = [
        n_hydroxy_pattern,
        n_di_hydroxy_pattern
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in n_hydroxy_patterns):
        return False, "No N-hydroxy modification found"

    return True, "Contains amino acid backbone with N-hydroxy modification"