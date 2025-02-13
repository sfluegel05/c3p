"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    A prenol is an alcohol possessing one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # An isoprene unit pattern represents the ongoing connectivity in prenol molecules
    many_isoprene_pattern = Chem.MolFromSmarts("C(C)(C)=C")
    
    if not mol.HasSubstructMatch(many_isoprene_pattern):
        return False, "No isoprene units found"

    # Check for terminal alcohol functionalities including OH, phosphates, and others
    term_alcohol_patterns = [
        Chem.MolFromSmarts("[OH]"),  # Direct OH group
        Chem.MolFromSmarts("COP(O)(=O)O"),  # Phosphate ester
        Chem.MolFromSmarts("C(=O)O"),  # Ester form often seen in derivatized forms
    ]
    
    for pattern in term_alcohol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains isoprene units with terminal alcohol group"

    return False, "No suitable terminal alcohol functionality found"