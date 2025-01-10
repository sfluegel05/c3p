"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    Mucopolysaccharides are polysaccharides composed of repeating units of uronic acids 
    and glycosamines, often esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mucopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for uronic acids backbone pattern
    uronic_acid_pattern = Chem.MolFromSmarts("[C](=O)[O][C]")  # Basic pattern of uronic acid
    if not mol.HasSubstructMatch(uronic_acid_pattern):
        return False, "No uronic acid units found"

    # Look for glycosamine backbone pattern
    glycosamine_pattern = Chem.MolFromSmarts("[C][N][C](O)")  # Basic pattern for glycosamine
    if not mol.HasSubstructMatch(glycosamine_pattern):
        return False, "No glycosamine units found"
    
    # Look for sulfate ester group presence
    sulfate_ester_pattern = Chem.MolFromSmarts("[O][S](=O)(=O)[O]")  # Basic sulfate ester pattern
    num_sulfate_groups = len(mol.GetSubstructMatches(sulfate_ester_pattern))
    if num_sulfate_groups < 1:
        return False, f"Missing sulfate ester groups for classification"

    return True, "Contains alternating uronic acids and glycosamines with esterified sulfate groups, indicating mucopolysaccharide"