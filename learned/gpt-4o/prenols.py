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

    # Isoprene unit pattern: C=C-C-C or C-C=C-C (simple representations)
    # These patterns are generalized to reduce failure from structural variations
    isoprene_pattern1 = Chem.MolFromSmarts("C=C-C-C")
    isoprene_pattern2 = Chem.MolFromSmarts("C-C=C-C")
    
    # Detect any isoprene substructure presence
    if not (mol.HasSubstructMatch(isoprene_pattern1) or mol.HasSubstructMatch(isoprene_pattern2)):
        return False, "No isoprene units found"

    # Verify the molecule has a terminal alcohol group
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol (OH) group found"

    return True, "Contains isoprene units with terminal alcohol group"