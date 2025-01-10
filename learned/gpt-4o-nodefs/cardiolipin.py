"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin molecule consists of a unique diphosphatidylglycerol core with four acyl chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for specific cardiolipin structural features
    # Pattern for the core cardiolipin structure: phosphatidylglycerol-phosphatidylglycerol
    # Use SMARTS patterns to identify known cardiolipin connectivity
    
    cardiolipin_core_pattern = Chem.MolFromSmarts("O(=P)(O)OCC(O)COP(O)(=O)O")  # Simplified core structure pattern
    
    if not mol.HasSubstructMatch(cardiolipin_core_pattern):
        return False, "Structure does not match cardiolipin core motif"

    # Check for exactly four ester linkages in correct positions signifying acyl chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Found {len(ester_matches)} ester linkages, need exactly 4 in cardiolipin-specific context"

    return True, "Contains structural motifs of cardiolipin (unique glycerol-phosphate core with four ester-linked acyl chains)"