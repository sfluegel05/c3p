"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is characterized by a nitrogen atom bonded to exactly 
    three carbon atoms, considering potential aromatic and aliphatic contexts.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Iterate through each atom in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check if the atom is nitrogen
            # List carbon atom neighbors
            carbon_neighbors = [neighbor for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
            
            # Check if the nitrogen has exactly three carbon neighbors
            if len(carbon_neighbors) == 3:
                # Ensure nitrogen has no other types of bonded atoms
                if all(neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
                    return True, "Molecule contains a tertiary amine (N bonded to 3 carbon atoms)"

    return False, "No tertiary amine found (N not bonded to 3 carbon atoms)"