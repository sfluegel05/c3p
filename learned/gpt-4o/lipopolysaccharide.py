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

    # Updated sugar patterns that are more likely in lipopolysaccharides.
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C(O)C1"),  # Glucose-like pyranose
        Chem.MolFromSmarts("O1COC(O)CC1"),  # General furanose
        Chem.MolFromSmarts("[C@H](O)[C@@H](O)[C@@H](O)C(O)C1OC(CO)C(O)1"),  # Specific repeating unit
    ]
    sugar_matches = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_matches:
        return False, "No appropriate sugar patterns detected"

    # Enhanced fatty acid patterns, considering three-hydroxytetradecanoic acid units.
    fatty_acid_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O"),  # Common long-chain acid
        Chem.MolFromSmarts("CCCCCCCCCCCCCCC(O)O"),  # Hydroxylated at terminal
        Chem.MolFromSmarts("CCC(O)CCCCCCCCCOC(=O)"),  # Similar to lipid A structures
    ]
    fatty_acid_matches = any(mol.HasSubstructMatch(pattern) for pattern in fatty_acid_patterns)
    if not fatty_acid_matches:
        return False, "Didn't match the specialized fatty acid chains like hydroxylated tetradecanoic"

    # Ensure presence of ester bonds, which are often vital in lipopolysaccharide backbones.
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 3:
        return False, f"Insufficient ester linkages; found {len(ester_matches)}"

    # Look for polysaccharide repeat units and phosphates which are part of lipopolysaccharide structure.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O[C]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Missing phosphate groups typical in lipopolysaccharides"

    return True, "Structure consistent with lipopolysaccharide characteristics"