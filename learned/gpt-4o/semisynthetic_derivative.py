"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Example natural product motifs such as steroid backbone, polyketide chains, etc.
    natural_product_patterns = [
        Chem.MolFromSmarts("C1CCC[C@]2(C1)CCC[C@H]3[C@H]2CC=C4C3=C(C=CC(=C4O)O)O"),  # Simplified steroid motif
        Chem.MolFromSmarts("C1=CC=CC=C1"),  # Aromatic ring (e.g., from polyketides)
    ]
    
    # Example of synthetic modifications like halogens, acylations, etc.
    synthetic_modifications_patterns = [
        Chem.MolFromSmarts("F"),  # Example: Fluorine, a common synthetically added group
        Chem.MolFromSmarts("Cl"), # Example: Chlorine, another common synthetic modification
        Chem.MolFromSmarts("C(=O)O"),  # ester function
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