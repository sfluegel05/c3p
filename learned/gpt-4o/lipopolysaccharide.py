"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    Lipopolysaccharides consist of a trisaccharide repeating unit with oligosaccharide side chains 
    and 3-hydroxytetradecanoic acid units.
    
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

    # Look for sugar patterns (repeating saccharide units)
    sugar_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C1")  # Basic pattern for sugars
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar patterns detected"

    # Look for fatty acid-like chain (specifically long chains like 3-hydroxytetradecanoic)
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(O)C(=O)O")  # heuristic pattern for tetradecanoic acid
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No long-chain fatty acid detected (e.g., 3-hydroxytetradecanoic acid)"

    # Check for multiple repeating units and branching typical of polysaccharides
    num_sugar_units = len(mol.GetSubstructMatches(sugar_pattern))
    if num_sugar_units < 3:
        return False, f"Insufficient sugar units for a lipopolysaccharide, found {num_sugar_units}"

    # As a simple implementation, further specificity checks like ensuring oligosaccharides or validating exact 
    # trisaccharide nature is beyond scope given complexity without detailed SMARTS for each feature
    # But assuming sugar and fatty acid presence is a fair first approximation for lipopolysaccharides

    return True, "Structure consistent with lipopolysaccharide characteristics"