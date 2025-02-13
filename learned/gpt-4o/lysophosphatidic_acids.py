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
    
    # Check for glycerol backbone with phosphate group attached
    glycerol_phosphate_pattern = Chem.MolFromSmarts("OCC(O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "Glycerol phosphate backbone not found"
    
    # Check for presence of at least one ester linkage attached to glycerol
    ester_pattern = Chem.MolFromSmarts("C(=O)OCC(O)CO")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No acyl chain esterified to the glycerol backbone"
    
    # Verify single ester linkage rule for simplicity; usually TAG structure
    # Look for presence of only one glycerol-attached carbon chain starting with C=O
    if len(ester_matches) > 1:
        return False, f"More than one ester linkage identified, found {len(ester_matches)}"

    # Additional checks can include typical lysophosphatidic acid characteristics like carbon chain lengths
    # All criteria met, it is a lysophosphatidic acid.
    return True, "The molecule is a lysophosphatidic acid"