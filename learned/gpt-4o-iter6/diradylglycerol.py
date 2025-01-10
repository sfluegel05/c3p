"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol has a glycerol backbone with exactly two substituent groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define glycerol backbones (C3O2 links accounted, not strictly linear to manage complexity)
    backbone_pattern = Chem.MolFromSmarts("[CX4](CO)(CO)O")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "Glycerol backbone pattern not matched"

    # Recognize ester (C=O) and ether (C-O) excluding rings more accurately
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[O]")
    ether_pattern = Chem.MolFromSmarts("[O][CX4]")

    # Recognize substituent links
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    # Only external substituent count is valid, immediate cycle neighbors should be omitted/included
    total_substituent_links = len(set(ester_matches + ether_matches))  # Ensures uniqueness and avoids duplicates

    # Ensure the molecular structure fits diradylglycerol profile (exactly 2 chains)
    if total_substituent_links != 2:
        return False, f"Expected exactly 2 substituent groups, found {total_substituent_links}"

    return True, "Valid diradylglycerol with two substituent groups connected via ester or ether bonds"