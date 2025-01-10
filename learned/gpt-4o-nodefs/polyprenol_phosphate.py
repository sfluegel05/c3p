"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is characterized by a polyprenyl chain terminated by a phosphate or diphosphate group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for isoprene units
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C-C")

    # Find matching substructures for isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    # Check for sufficient isoprene units (polyprenol usually has multiple units)
    if len(isoprene_matches) < 2: # Adjust this number appropriately to classify based on definition
        return False, "Too few isoprene units for polyprenol"

    # Define SMARTS patterns for phosphate and diphosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[O]P(=O)(O)O")
    diphosphate_pattern = Chem.MolFromSmarts("OP(=O)([O])OP(=O)([O])O")

    # Check whether the molecule contains a terminal phosphate or diphosphate group
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)

    if not (phosphate_matches or diphosphate_matches):
        return False, "No terminal phosphate or diphosphate group found"

    return True, "Contains polyprenyl chain with terminal phosphate or diphosphate group"