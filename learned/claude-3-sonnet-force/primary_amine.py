"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:33842 primary amine
A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get nitrogen atoms
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

    # Check if there is at least one nitrogen atom
    if not n_atoms:
        return False, "No nitrogen atom found"

    # Check for primary amine groups
    primary_amine_found = False
    for n_atom in n_atoms:
        # Check if nitrogen has exactly one hydrogen and one carbon neighbor
        h_neighbors = sum(1 for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 1)
        c_neighbors = sum(1 for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
        if h_neighbors == 1 and c_neighbors == 1:
            # This is a primary amine nitrogen
            primary_amine_found = True
            break

    if not primary_amine_found:
        return False, "No primary amine group found"

    return True, "Contains a primary amine group (-NH2) bonded to a carbon"