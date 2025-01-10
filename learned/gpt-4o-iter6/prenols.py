"""
Classifies: CHEBI:26244 prenols
"""
from rdkit import Chem

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with the formula H-[CH2C(Me)=CHCH2]nOH, 
    containing one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined pattern to capture isoprene units with possible geometric isomerism
    isoprene_smarts = "C(=C(C)C)[CH2]"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)
    
    # Search for at least one isoprene unit
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 1:
        return False, "No isoprene units found"
    
    # Check for terminal alcohol - more specific to ensure the alcohol is terminal
    # A terminal alcohol will have one OH group bonded to a non-carbonyl carbon and end of chain
    terminal_alcohol_pattern = Chem.MolFromSmarts("[C](O)")
    terminal_alcohol_matches = mol.GetSubstructMatches(terminal_alcohol_pattern)
    
    # We need to be 100% sure that this is indeed terminal
    is_terminal = False
    for match in terminal_alcohol_matches:
        alcohol_atom = mol.GetAtomWithIdx(match[1])
        # Ensure that this OH is terminal
        if alcohol_atom.GetDegree() == 1:
            is_terminal = True
            break
    
    if not is_terminal:
        return False, "Missing or non-terminal alcohol group"
    
    # Basic logic now checks for presence of isoprene unit(s) and terminal OH group
    return True, f"Contains {len(isoprene_matches)} isoprene units with a terminal alcohol group"

__metadata__ = {'chemical_class': {'name': 'prenol'}}