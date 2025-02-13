"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:51379 polyamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is an organic amino compound that contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and standardize molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    
    # Count amino groups
    amino_pattern = Chem.MolFromSmarts("[NH2]")
    amino_groups = len(mol.GetSubstructMatches(amino_pattern))
    
    # Polyamines must have at least 2 amino groups
    if amino_groups < 2:
        return False, "Fewer than 2 amino groups"
    
    # Check for specific substructures/patterns characteristic of polyamines
    polyamine_patterns = [
        Chem.MolFromSmarts("[NH2]-[CH2]-[CH2]-[CH2]-[NH2]"),  # Alkyl chain with amino groups
        Chem.MolFromSmarts("[NH2]-[c]-[c]-[c]-[NH2]"),  # Aromatic ring with amino groups
        Chem.MolFromSmarts("[NH2]-[c]-[c]-[NH2]"),  # Aromatic ring with amino groups
        Chem.MolFromSmarts("[NH2]-[CH2]-[CH2]-[NH2]"),  # Ethylenediamine
        Chem.MolFromSmarts("[NH2]-[CH2]-[NH]-[CH2]-[NH2]"),  # Triazine with amino groups
        # Add more patterns as needed
    ]
    
    for pattern in polyamine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {amino_groups} amino groups and matches polyamine substructure pattern"
    
    # Check for specific exclusion patterns
    exclusion_patterns = [
        Chem.MolFromSmarts("[N+]"),  # Quaternary amines
        Chem.MolFromSmarts("[N-]"),  # Anionic amines
        Chem.MolFromSmarts("[N+]=[N-]"),  # Azo compounds
        Chem.MolFromSmarts("[N-]=[N+]=[N-]"),  # Azide compounds
        Chem.MolFromSmarts("[N]=[N]"),  # Azo compounds
        Chem.MolFromSmarts("[N]#[C]"),  # Nitriles
        Chem.MolFromSmarts("[N]=[C]"),  # Amides
        # Add more patterns as needed
    ]
    
    for pattern in exclusion_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, f"Contains {amino_groups} amino groups but matches exclusion pattern"
    
    return True, f"Contains {amino_groups} amino groups and does not match any exclusion pattern"