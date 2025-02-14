"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is derived by replacing one hydrogen atom of ammonia with a hydrocarbyl group.

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
    
    # Loop through all nitrogen atoms to check for primary amine characteristics
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            num_hydrogens = atom.GetTotalNumHs()
            num_neighbors = len(atom.GetNeighbors())

            # Check if this nitrogen forms a primary amine
            if num_hydrogens == 2 and num_neighbors == 1:
                neighbor = atom.GetNeighbors()[0]
                
                # The nitrogen must be connected to a carbon atom
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    return True, "Primary amine group (NH2) attached to carbon detected"

    return False, "No primary amine group (NH2) attached to carbon detected"