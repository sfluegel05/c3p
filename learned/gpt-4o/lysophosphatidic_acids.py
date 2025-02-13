"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
from rdkit import Chem

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.

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
    
    # Check for glycerol backbone with one ester linked fatty acid
    glycerol_pattern = Chem.MolFromSmarts("OC[C@H](O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"
    
    # Check for an ester linkage (C(=O)OC) indicating a single acyl chain
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected one ester linkage, found {len(ester_matches)}"
    
    # Check for phosphate group present in lysophosphatidic acids
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found"

    # If all checks passed
    return True, "The molecule is a lysophosphatidic acid"