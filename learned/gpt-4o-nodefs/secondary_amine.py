"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine has a nitrogen atom bonded to exactly two carbon atoms.

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

    # Iterate over all atoms to find secondary amine characteristic
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check for Nitrogen atom
            neighbors = atom.GetNeighbors()
            num_carbon = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 6)
            num_nonH = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() != 1)

            # A secondary amine has Nitrogen connected to two Carbon atoms and an optional hydrogen or other atom
            if num_carbon == 2 and num_nonH == 2:
                return True, "Contains a secondary amine group"
    
    return False, "Does not contain a secondary amine group"

# You can test the function with a SMILES string to see if it correctly classifies secondary amines.