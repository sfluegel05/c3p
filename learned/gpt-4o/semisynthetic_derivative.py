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
        Chem.MolFromSmarts("C1CCC[C@]2(C1)CCC[C@H]3[C@H]2CC=C4C3=C(C=CC(=C4O)O)O"),  # Steroid-like
        Chem.MolFromSmarts("[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@H](C1)O"),  # Sugar-like
        Chem.MolFromSmarts("C=C[C@H](C)[C@H](O)[C@@H](C)C=C"),  # Terpene-like structure
        Chem.MolFromSmarts("CN1C=NC2=C1C(=O)N(C(=O)N2)C"),  # Nucleoside base-like motif
        Chem.MolFromSmarts("[C@H]1C=C[C@H](O)[C@@H](O)C1"),  # Aromatic rings with hydroxyls
    ]
    
    # Expanded set of synthetic modifications
    synthetic_modifications_patterns = [
        Chem.MolFromSmarts("F"),  # Fluorine
        Chem.MolFromSmarts("[Cl]"), # Chlorine
        Chem.MolFromSmarts("C(=O)O"),  # Ester function
        Chem.MolFromSmarts("N(C)(C)"), # Alkylated amine
        Chem.MolFromSmarts("C#N"),  # Nitrile group
        Chem.MolFromSmarts("S(=O)(=O)"),  # Sulfonyl group
        Chem.MolFromSmarts("[C&R2]-[C&!R]"),  # Non-ring alkyl chain
        Chem.MolFromSmarts("[NH2]-C(=O)"),  # Amide linkages
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