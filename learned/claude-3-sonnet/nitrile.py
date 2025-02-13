"""
Classifies: CHEBI:18379 nitrile
"""
"""
Classifies: CHEBI:46419 nitrile
A compound having the structure RC#N; thus a C-substituted derivative of hydrocyanic acid, HC#N. 
In systematic nomenclature, the suffix nitrile denotes the triply bound #N atom, not the carbon atom attached to it.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nitrile(smiles: str):
    """
    Determines if a molecule is a nitrile based on its SMILES string.
    A nitrile contains a cyano group (-C#N) attached to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyano group (-C#N) anywhere in the molecule
    cyano_pattern = Chem.MolFromSmarts("[C#N]")
    if not mol.HasSubstructMatch(cyano_pattern):
        return False, "No cyano group (-C#N) found"
    
    # Check for common non-nitrile functional groups
    non_nitrile_patterns = [
        Chem.MolFromSmarts("[C](=O)[N]"),  # Amide
        Chem.MolFromSmarts("[C](=O)[O]"),  # Ester
        Chem.MolFromSmarts("[N+](=O)[O-]"),  # Nitro
        Chem.MolFromSmarts("[N+]#[C-]"),  # Isocyanide
        Chem.MolFromSmarts("[N,O][S](=O)=O"),  # Sulfonamide, sulfonate
        Chem.MolFromSmarts("[N]C(=O)C"),  # Amide
        Chem.MolFromSmarts("[N]C(=O)N"),  # Urea
        Chem.MolFromSmarts("[N]C(=O)O")  # Carbamate
    ]
    
    for pattern in non_nitrile_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, f"Contains a {pattern.GetSmarts()} group, which is not a nitrile"
    
    # Additional checks for common nitrile structural features
    if mol.GetNumAtoms() < 3:
        return False, "Too small to be a nitrile compound"
    
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7) > 2:
        return False, "Too many nitrogen atoms for a typical nitrile"
    
    # If cyano group present and no common non-nitrile groups, classify as nitrile
    return True, "Contains a cyano group (-C#N) attached to a carbon atom"