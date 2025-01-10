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
    
    # Define a more complex isoprene unit pattern to account for variation and repetition
    isoprene_pattern = Chem.MolFromSmarts("[C;R0](=C)[C;R0]-C-C")
    
    # Define more flexible terminal alcohol pattern considering possible connectivity
    terminal_alcohol_pattern = Chem.MolFromSmarts("[C;R0][OH]")
    
    # Search for isoprene units
    isoprene_blocks = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_blocks) < 1:
        return False, "No or insufficient isoprene units found"
    
    # Check alcohol presence and ensure it is terminal
    alcohol_matches = mol.GetSubstructMatches(terminal_alcohol_pattern)
    if len(alcohol_matches) == 0:
        return False, "No or non-terminal alcohol group found"

    # Validate the position of OH to be at the terminus of the molecule
    terminal_endings = [match for match in alcohol_matches if mol.GetAtomWithIdx(match[0]).GetDegree() == 1]
    if len(terminal_endings) == 0:
        return False, "Alcohol group not at the terminal end"

    return True, "Contains isoprene units with a terminal alcohol group"

__metadata__ = {'chemical_class': {'name': 'prenol'}}