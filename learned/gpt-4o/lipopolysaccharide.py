"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides are complex molecules including oligosaccharides, fatty acids, and are a major
    component of Gram-negative bacteria cell walls.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Additional sugar patterns for more comprehensive detection
    sugar_patterns = [
        Chem.MolFromSmarts("O1[C@@H]([C@@H](O)[C@H](O)[C@H]([C@@H]1O)CO)C(=O)O"),  # Generic hexose-like structures
        Chem.MolFromSmarts("C1=CC(=O)OC(=C1)O"),  # Aromatic hexose
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),  # Simplified pyranose
    ]
    sugar_matches = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_matches:
        return False, "No known lipopolysaccharide sugar patterns detected"

    # Refine fatty acid chain detection for hydroxytetradecanoic structures
    fatty_acid_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O"),  # Long-chain acids
        Chem.MolFromSmarts("CCC(O)CCCCCCCCCC(=O)OC"),  # Hydroxy-acyl variants
    ]
    fatty_acid_matches = any(mol.HasSubstructMatch(pattern) for pattern in fatty_acid_patterns)
    if not fatty_acid_matches:
        return False, "No recognized fatty acid chain patterns"

    # Look for ester and glycosidic linkages as indicative of lipopolysaccharide backbone
    glycosidic_pattern = Chem.MolFromSmarts("OC1CC(O)C(O)C(O)C1") 
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Missing core glycosidic linkages"

    # Check for phosphate groups integral to lipopolysaccharide structure
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups detected"

    return True, "Structural features consistent with lipopolysaccharide"