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

    # Identify a generalized glycerol-phosphate linkage pattern allowing stereochemistry variations
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)([O-,O])O")

    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol-phosphate linkage found"
        
    # Enhanced pattern for acyl chain with variability in chain and unsaturation
    acyl_chain_pattern = Chem.MolFromSmarts("C(=O)OC")
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)

    # Ensure exactly one acyl group matches required ester linkages
    if len(acyl_matches) != 1:
        return False, f"Found {len(acyl_matches)} acyl chains, need exactly 1"

    # Additional checks can be added based on false positives and specific exceptions
    # Ensure structure doesn't match exclusions (complex lipids or non-related phosphates)
    
    return True, "Matches all structural features of a lysophosphatidic acid"