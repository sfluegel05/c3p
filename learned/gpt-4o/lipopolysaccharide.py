"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Look for more general sugar patterns including various linkages (five or six-membered rings)
    sugar_patterns = [
        Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1"),  # Basic glucose
        Chem.MolFromSmarts("O1COC(O)C1"),  # Furanose
        Chem.MolFromSmarts("O1C(COC(O)C1)O"),  # Pyranose
    ]
    sugar_matches = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not sugar_matches:
        return False, "No sugar patterns detected"

    # Look for long-chain fatty acid patterns (explore different chain lengths and hydroxylation)
    fatty_acid_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O"),  # simple long acid
        Chem.MolFromSmarts("CCCCCCCCCCCCC(O)C(=O)O"),  # hydroxylated fatty acid
        Chem.MolFromSmarts("CCC(O)CCCCCCCCCC=O"),  # alternative hydroxyl group
    ]
    fatty_acid_matches = any(mol.HasSubstructMatch(pattern) for pattern in fatty_acid_patterns)
    if not fatty_acid_matches:
        return False, "No long-chain fatty acid detected, e.g., hydroxylated fatty acids"

    # Check for evidence of polysaccharide structure, ensuring branching or repeating nature
    num_sugar_units = sum(len(mol.GetSubstructMatches(pattern)) for pattern in sugar_patterns)
    if num_sugar_units < 3:
        return False, f"Insufficient sugar units for a lipopolysaccharide, found {num_sugar_units}"

    # Look for ester linkages typical of lipopolysaccharide lipid A groups
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing ester linkages typical for lipid structures"

    return True, "Structure consistent with lipopolysaccharide characteristics"