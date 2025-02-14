"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine has a nitrogen atom connected to two carbon atoms and one hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Atomic number of nitrogen
            # Get the neighbor atoms of the nitrogen
            neighbors = atom.GetNeighbors()
            carbon_count = 0
            hydrogen_count = 0

            # Count carbon and hydrogen neighbors
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    carbon_count += 1
                elif neighbor.GetAtomicNum() == 1:  # Hydrogen
                    hydrogen_count += 1

            # Check if nitrogen has exactly two carbon and one hydrogen neighbors
            if carbon_count == 2 and hydrogen_count == 1:
                return True, "Nitrogen is bonded to two carbons and one hydrogen, indicative of a secondary amine"

    return False, "No nitrogen atom bonded to exactly two carbons and one hydrogen found"