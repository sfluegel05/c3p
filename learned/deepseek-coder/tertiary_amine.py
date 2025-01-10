"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
A tertiary amine is a compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three carbon atoms (hydrocarbyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Check if the nitrogen is bonded to exactly 3 carbon atoms
            carbon_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
            if len(carbon_neighbors) == 3:
                return True, "Contains a nitrogen atom bonded to three carbon atoms (tertiary amine)"

    # If no such nitrogen is found, it's not a tertiary amine
    return False, "No nitrogen atom bonded to three carbon atoms found"

# Example usage:
# print(is_tertiary_amine("CCN(CC)CC"))  # Should return (True, "Contains a nitrogen atom bonded to three carbon atoms (tertiary amine)")