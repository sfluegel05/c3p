"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compounds (compounds with at least one carbon-chlorine bond)
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound by checking for any carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has at least one C-Cl bond (any bond type), False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # Check if atom is chlorine
            # Check all neighboring atoms to see if any are carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon has atomic number 6
                    return True, "Contains at least one carbon-chlorine bond"
    
    # If no C-Cl bonds found
    return False, "No carbon-chlorine bonds found"