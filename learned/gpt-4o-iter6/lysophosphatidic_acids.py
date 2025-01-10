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

    # Identify a more general glycerol-phosphate linkage pattern without strict stereochemistry
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)O")
    
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol-phosphate linkage found"
        
    # Use a more specific pattern refined for a single acyl chain
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)OC")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    
    # Accounting for ester linkages, ensure exactly one acyl group is present
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl chains, need exactly 1"

    # Further checks can be implemented here for specificity

    # If all conditions are met, it is classified as a lysophosphatidic acid
    return True, "Matches all structural features of a lysophosphatidic acid"