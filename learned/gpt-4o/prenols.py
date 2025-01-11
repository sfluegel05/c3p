"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol contains one or more isoprene units and an alcohol (OH) group
    or a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification or exclusion
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined isoprene unit patterns
    # Common isoprene configurations - linear or terminal methyl as seen in prenol forms
    isoprene_unit_1 = Chem.MolFromSmarts("C(=C)C(C)C")  # Example pattern for typical isoprene
    isoprene_unit_2 = Chem.MolFromSmarts("CC(=CC)C")    # Alternate pattern
    
    # Searching both isoprene patterns 
    num_isoprene_units = len(mol.GetSubstructMatches(isoprene_unit_1)) + len(mol.GetSubstructMatches(isoprene_unit_2))

    # Check for isoprene units
    if num_isoprene_units == 0:
        return False, "No isoprene units found"

    # Check for alcohol group
    oh_group = Chem.MolFromSmarts("[OX2H]")
    has_oh_group = mol.HasSubstructMatch(oh_group)
    
    # Check for phosphate groups
    phosphate_group = Chem.MolFromSmarts("P(=O)(O)[O-]")    # Refine to reflect divalent portrayal and deoxidized forms
    has_phosphate_group = mol.HasSubstructMatch(phosphate_group)

    # Ensure at least one of the groups (OH/phosphate) is present
    if has_oh_group:
        return True, f"Contains {num_isoprene_units} isoprene unit(s) and an alcohol group (OH)"
    elif has_phosphate_group:
        return True, f"Contains {num_isoprene_units} isoprene unit(s) and a phosphate group"
    
    return False, "No alcohol (OH) or phosphate group found"