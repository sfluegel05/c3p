"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS pattern: Acyl group attached at C1 with intact secondary alcohols
    # We identify C1-O-CO-[C] structure and ensure it's exclusive to this position
    # with a secondary alcohol configuration remaining nearby

    # Pattern for ester at first carbon including stereochemistry/sn-glycerol form:
    monoglyceride_pattern = Chem.MolFromSmarts('OCC(OC=O)[C@@H](O)CO')  # Matching glycerol ester at primary position

    # Check for presence of the specific pattern
    if not mol.HasSubstructMatch(monoglyceride_pattern):
        return False, "1-monoglyceride structural pattern not matched"

    # Verify esterification is at the primary position of glycerol's backbone
    ester_pattern = Chem.MolFromSmarts('O=C[O][CH2][CH](O)CO')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Acyl linkage not specifically at position 1"
    
    # Additional: Ensure there's no other ester linkage taking away features
    branch_matches = mol.GetSubstructMatches(monoglyceride_pattern)
    if len(branch_matches) > 1:
        return False, f"Multiple esterification detected, expected only one at the primary position"

    return True, "Verified 1-monoglyceride with correct acyl linkage at the position 1"