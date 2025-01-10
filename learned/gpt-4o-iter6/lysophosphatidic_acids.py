"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid is a monoacylglycerol phosphate with specific structural features.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify a more general glycerol-phosphate linkage pattern
    # Variations include different stereochemistry
    glycerol_phosphate_patterns = [
        Chem.MolFromSmarts("OCC(O)COP(=O)(O)O"),  # without stereochemistry
        Chem.MolFromSmarts("O[C@H](CO)COP(=O)(O)O"),  # specific stereochemistry
        Chem.MolFromSmarts("O[C@@H](CO)COP(=O)(O)O")   # specific opposite stereochemistry
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in glycerol_phosphate_patterns):
        return False, "No glycerol-phosphate linkage found"
        
    # Use a broad pattern to check for exactly one acyl chain
    # Allow for chain variation and different carbon lengths
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)[CX4,CX3]")  # Simple ester linkage 
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl chains, need exactly 1"

    # If all conditions are met, it is classified as a lysophosphatidic acid
    return True, "Matches all structural features of a lysophosphatidic acid"