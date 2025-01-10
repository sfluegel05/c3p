"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar contains a carbohydrate backbone with one or more oxygens
    or hydroxy groups replaced by sulfur or -SR groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbohydrate pattern: 5 or 6 membered oxygen-containing ring
    carbohydrate_pattern = Chem.MolFromSmarts('[C@H1]1([C@@H1]([C1@H]O)O)O')
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No carbohydrate-like structure found"
        
    # Look for sulfur substitution patterns in place of oxygen
    sulfur_pattern_1 = Chem.MolFromSmarts('[C@H1]1([C@@H1]([C1@H]S)O)O')  # Sulfur in ring
    sulfur_pattern_2 = Chem.MolFromSmarts('SC')  # Simple -SR substitution pattern

    if mol.HasSubstructMatch(sulfur_pattern_1) or mol.HasSubstructMatch(sulfur_pattern_2):
        return True, "Sulfur substitution found in carbohydrate-like structure"

    return False, "No sulfur substitution in expected carbohydrate structure locations"