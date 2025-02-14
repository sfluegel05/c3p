"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is likely to be a mucopolysaccharide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has characteristics of a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enhanced check for uronic acid-like features
    uronic_acid_patterns = [
        Chem.MolFromSmarts("[C;R0](=O)[O;R0]"),  # Carboxyl group, often in uronic acids
        Chem.MolFromSmarts("[C;R0](=O)O[CH1][CH2]")  # Sugars leading to uronic acids often have this attachment
    ]
    
    # Check for some nitrogen functionalities hinting glycosamines (fairly generalized)
    glycosamine_patterns = [
        Chem.MolFromSmarts("[*][NX3][CX4]"),
        Chem.MolFromSmarts("[*][NX3][CX3]=O")
    ]
    
    # Check for sulfate esterification
    sulfate_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])([OX2])")

    # Checking for individual and alternating patterns
    has_uronic_acid = any(mol.HasSubstructMatch(pat) for pat in uronic_acid_patterns)
    has_glycosamine = any(mol.HasSubstructMatch(pat) for pat in glycosamine_patterns)
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)

    if has_uronic_acid and has_glycosamine:
        if has_sulfate:
            return True, "Characteristics of mucopolysaccharide, including sulfate esterification"
        else:
            return True, "Characteristics of mucopolysaccharide without detectable sulfate groups"
        
    return False, "Missing key features of a mucopolysaccharide"