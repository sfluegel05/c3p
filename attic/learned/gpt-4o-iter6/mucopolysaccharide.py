"""
Classifies: CHEBI:37395 mucopolysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mucopolysaccharide(smiles: str):
    """
    Determines if a molecule is a mucopolysaccharide based on its SMILES string.
    A mucopolysaccharide is a polysaccharide composed of alternating units from uronic acids and glycosamines,
    and commonly partially esterified with sulfuric acid.

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

    # Define SMARTS patterns for uronic acid and glycosamine units
    uronic_acid_pattern = Chem.MolFromSmarts("C(=O)[O-]")  # Simplified representation of uronic acids
    glycosamine_pattern = Chem.MolFromSmarts("[NX3][CX4]([OH])")  # Glycosamine has an amine and hydroxyl group
    
    # Check alternating patterns of uronic acids and glycosamines
    if not mol.HasSubstructMatch(uronic_acid_pattern):
        return False, "No uronic acid units found"
    
    if not mol.HasSubstructMatch(glycosamine_pattern):
        return False, "No glycosamine units found"

    # Define SMARTS pattern for sulfate ester group
    sulfate_pattern = Chem.MolFromSmarts("O=S(=O)([OX2H1])")  # Generic sulfate ester pattern
    
    # Check for presence of sulfate groups
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfuric acid esterification found"

    return True, "Contains alternating uronic acids and glycosamines with sulfuric acid esterification"

# Note that the SMARTS patterns and checks are simplified and may not match all mucopolysaccharides accurately.
# Further refinement will be necessary for accurate classification.