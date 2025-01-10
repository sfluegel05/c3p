"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with the formula H-[CH2C(Me)=CHCH2]nOH, containing one or more isoprene units.

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
    
    # Pattern for isoprene unit
    isoprene_pattern = Chem.MolFromSmarts("C(C)=C(CCC)=C")
    
    # Search for isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No or insufficient isoprene units found"
    
    # SMARTS for terminal alcohol group
    terminal_alcohol_pattern = Chem.MolFromSmarts("[C;!$(C=O)][OH]")
    
    # Check for a terminal alcohol group
    if not mol.HasSubstructMatch(terminal_alcohol_pattern):
        return False, "Missing or non-terminal alcohol group"
    
    # Assuming isoprene and OH presence means prenol
    return True, "Contains isoprene units with terminal alcohol group"

__metadata__ = {'chemical_class': {'name': 'prenol'}}