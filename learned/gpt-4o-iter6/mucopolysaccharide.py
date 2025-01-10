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

    # Refined SMARTS pattern for uronic acid (R-C(=O)[O-] where R is a carbohydrate portion)
    uronic_acid_pattern = Chem.MolFromSmarts("OC(=O)C[OH]")  # R-O-C(=O)C with a hydroxyl on the sugar

    # Refined SMARTS pattern for glycosamine (sugar unit with amine)
    glycosamine_pattern = Chem.MolFromSmarts("C(C([OH])C[OH])[NX3]")  # Simplified sugar-like pattern with amine
    
    # Define SMARTS pattern for sulfate ester group
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O[CX4]")  # Sulfate ester pattern bonded to carbon

    # Check for alternating uronic acid and glycosamine units
    uronic_acid_matches = mol.GetSubstructMatches(uronic_acid_pattern)
    glycosamine_matches = mol.GetSubstructMatches(glycosamine_pattern)

    if not uronic_acid_matches or not glycosamine_matches:
        return False, "No alternating sequence of uronic acids and glycosamines found"

    # Check for sulfate esterification
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfuric acid esterification found"

    # Additional checks for alternation pattern could involve checking adjacency of these matches
    uronic_positions = {match[0] for match in uronic_acid_matches}
    glyco_positions = {match[0] for match in glycosamine_matches}
    
    # Check adjacency or direct bond connections (conceptual demonstration, not exact check here)
    alternating = False
    for ua_pos in uronic_positions:
        for ga_pos in glyco_positions:
            if mol.GetBondBetweenAtoms(ua_pos, ga_pos):
                alternating = True

    if not alternating:
        return False, "No alternating sequence of uronic acids and glycosamines found"

    return True, "Contains alternating uronic acids and glycosamines with sulfuric acid esterification"