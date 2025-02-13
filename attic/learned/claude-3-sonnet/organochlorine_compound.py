"""
Classifies: CHEBI:36683 organochlorine compound
"""
"""
Classifies: organochlorine compound
Definition: A compound containing at least one carbon-chlorine bond
"""
from rdkit import Chem

def is_organochlorine_compound(smiles: str):
    """
    Determines if a molecule is an organochlorine compound based on its SMILES string.
    An organochlorine compound must have at least one carbon-chlorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organochlorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find chlorine atoms
    chlorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17]
    
    if not chlorine_atoms:
        return False, "No chlorine atoms found"

    # Check if any chlorine is bonded to carbon
    for cl_atom in chlorine_atoms:
        for neighbor in cl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atomic number
                return True, "Contains carbon-chlorine bond"
                
    return False, "Contains chlorine but no carbon-chlorine bonds"