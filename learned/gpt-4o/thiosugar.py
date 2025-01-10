"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar has a carbohydrate backbone with one or more oxygens or hydroxy
    groups replaced by sulfur or -SR groups.

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
    
    # More flexible pattern for carbohydrate recognition 
    # Allow for ring or open-chain forms with multiple hydroxyls
    carbohydrate_pattern = Chem.MolFromSmarts("C(O)(C(=O)O)O")
    
    # Check if there is at least one recognizable carbohydrate motif
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No recognizable carbohydrate motifs found"
    
    # Look for sulfur substituents replacing oxygens or hydroxyls
    sulfur_substitution_patterns = [
        Chem.MolFromSmarts("[SX2]"),       # Sulfide
        Chem.MolFromSmarts("[SX3]=O"),    # Sulfoxide
        Chem.MolFromSmarts("[SX4](=O)(=O)"), # Sulfone
        Chem.MolFromSmarts("S-[C]")          # Sulfur linked to carbon indicating an -SR group
    ]
    
    for pattern in sulfur_substitution_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Sulfur substitution found in potential carbohydrate structure"
    
    return False, "No sulfur substitution found in expected carbohydrate structure locations"