"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define SMARTS for polyprenyl chain (patterns of isoprene units)
    isoprene_pattern = Chem.MolFromSmarts("C=C-C-C")
    
    # Find matching substructures for isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:  # Must be more than one isoprene unit
        return False, "Too few isoprene units for polyprenol"

    # Define SMARTS for phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    diphosphate_pattern = Chem.MolFromSmarts("P(=O)(O)OP(=O)(O)O")
    
    # Check if it contains a phosphate or diphosphate group
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    has_diphosphate = mol.HasSubstructMatch(diphosphate_pattern)
    
    if not (has_phosphate or has_diphosphate):
        return False, "No phosphate or diphosphate group found"

    return True, "Contains polyprenyl chain with terminal phosphate or diphosphate group"