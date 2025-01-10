"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem

def is_semisynthetic_derivative(smiles: str):
    """
    Determines if a molecule is a semisynthetic derivative based on its SMILES string.
    A semisynthetic derivative is an organic molecular entity derived from a natural product
    by partial chemical synthesis.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a semisynthetic derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded set of natural product motifs
    natural_product_patterns = [
        Chem.MolFromSmarts("[C@H]1[C@@H]([C@H]([C@@H](C1)O)O)O"),  # Simple sugar/furanose-like
        Chem.MolFromSmarts("C1CCC[C@H]1"),  # Cyclohexane-like - representing a more generalized motif
        Chem.MolFromSmarts("C=C[C@H](C)C"),  # Simple terpene-like backbone
    ]
    
    # Expanded set of synthetic modifications
    synthetic_modifications_patterns = [
        Chem.MolFromSmarts("F"),  # Fluorine
        Chem.MolFromSmarts("[Cl]"), # Chlorine, presence in drug-like or synthetic contexts
        Chem.MolFromSmarts("C(=O)O"),  # Ester function more specifically within certain contexts
        Chem.MolFromSmarts("N(C)C(=O)"), # Amides with specific alkylation
        Chem.MolFromSmarts("/C=C\\C(=O)"),  # Vinyl ketones; potential semisynthetic modification
    ]

    # Check for natural product motifs
    has_natural_motif = any(mol.HasSubstructMatch(pattern) for pattern in natural_product_patterns)
    
    # Check for synthetic modifications
    has_synthetic_modification = any(mol.HasSubstructMatch(pattern) for pattern in synthetic_modifications_patterns)
    
    # Determine if both conditions are met
    if has_natural_motif and has_synthetic_modification:
        return True, "Contains both natural product motifs and synthetic modifications"
    elif not has_natural_motif:
        return False, "No recognizable natural product motif found"
    elif not has_synthetic_modification:
        return False, "No recognizable synthetic modification found"
    
    return False, "Does not meet criteria for semisynthetic derivative"