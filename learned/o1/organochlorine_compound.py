"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound contains at least one carbon-chlorine (C-Cl) bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all chlorine atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 17:  # Atomic number of chlorine
            # Check if any neighbor atom is carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Atomic number of carbon
                    return True, "Contains at least one carbon-chlorine (C-Cl) bond"

    # If no C-Cl bonds are found
    return False, "No carbon-chlorine (C-Cl) bonds found"