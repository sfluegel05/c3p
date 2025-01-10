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
    
    # Improved pattern for core cardiolipin structure
    cardiolipin_core_pattern = Chem.MolFromSmarts("O=P(O)(OCC(O)COP(O)(O)CC(O)CO)O")
    
    if not mol.HasSubstructMatch(cardiolipin_core_pattern):
        return False, "Structure does not match cardiolipin core motif"

    # Check for exactly four ester linkages representing the acyl chains
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]CO")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Found {len(ester_matches)} ester linkages, need at least 4 in cardiolipin-specific context"
    
    # Ensure phosphates link through glycerol (two per cardiolipin)
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)OCC(O)CO")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 2:
        return False, f"Found {len(phosphate_matches)} phosphate linkages, need exactly 2 for cardiolipin"

    return True, "Contains structural motifs of cardiolipin (unique glycerol-phosphate core with four ester-linked acyl chains)"