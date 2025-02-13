"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is derived by replacing one hydrogen atom of ammonia with a hydrocarbyl group,
    typically resulting in an "NH2" group attached to a carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a primary amine, False otherwise
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
            carbon_attached = 0

            # Check if this nitrogen forms a primary amine
            if num_hydrogens == 2:
                # Make sure the nitrogen is attached to exactly one carbon (primary amine condition)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Carbon
                        carbon_attached += 1

            if carbon_attached == 1 and num_hydrogens == 2 and num_neighbors == 1:
                return True, "Primary amine group (NH2) attached to carbon detected"

    return False, "No primary amine group (NH2) attached to carbon detected"