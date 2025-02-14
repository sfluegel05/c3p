"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is a compound containing at least one carbon atom and at least one fluorine atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Initialize flags for carbon and fluorine presence
    has_carbon = False
    has_fluorine = False
    
    # Iterate over all atoms to check for carbon and fluorine
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 6:
            has_carbon = True
        elif atomic_num == 9:
            has_fluorine = True
        # If both are found, no need to continue
        if has_carbon and has_fluorine:
            return True, "Contains both carbon and fluorine atoms"
    
    # Check which atom was missing
    if not has_carbon and not has_fluorine:
        return False, "Contains neither carbon nor fluorine atoms"
    elif not has_carbon:
        return False, "No carbon atoms found"
    elif not has_fluorine:
        return False, "No fluorine atoms found"
    
    return False, "Does not meet criteria for an organofluorine compound"