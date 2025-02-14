"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:35322 organofluorine compound
An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carbon-fluorine bonds using SMARTS pattern
    c_f_pattern = Chem.MolFromSmarts("[C-F]")
    matches = mol.GetSubstructMatches(c_f_pattern)
    
    # If any carbon-fluorine bond is found, classify as organofluorine compound
    if matches:
        return True, "Contains at least one carbon-fluorine bond"
    
    # Check for special cases: trifluoromethyl and other fluorine-containing groups not bonded to carbon
    excluded_patterns = [
        Chem.MolFromSmarts("[C-N-F]"),      # Fluorine bonded to nitrogen
        Chem.MolFromSmarts("[C-O-F]"),      # Fluorine bonded to oxygen
        Chem.MolFromSmarts("[$(C(F)(F)F)]") # Trifluoromethyl group
    ]
    
    for pattern in excluded_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains non-carbon-fluorine bonds or functional groups"
    
    # If no carbon-fluorine bond is found, classify as non-organofluorine compound
    return False, "No carbon-fluorine bond found"