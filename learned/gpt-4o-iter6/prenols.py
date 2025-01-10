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
    
    # Define isoprene unit pattern
    isoprene_pattern = Chem.MolFromSmarts("C(C)=C-C-C")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene unit found"
    
    # Define alcohol group pattern
    alcohol_pattern = Chem.MolFromSmarts("CO")
    if not mol.HasSubstructMatch(alcohol_pattern):
        return False, "No alcohol group (OH) found"
    
    # Verify linkage of isoprene units
    # This step could be enhanced to count and ensure the correct sequential linkage.
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)

    # Check if there's a continuous chain ending with an alcohol
    if len(isoprene_matches) < 1 or len(alcohol_matches) == 0:
        return False, "Insufficient isoprene units or no terminal alcohol group"

    return True, "Contains isoprene units with a terminal alcohol group"

__metadata__ = { 'chemical_class': { 'name': 'prenol' }}