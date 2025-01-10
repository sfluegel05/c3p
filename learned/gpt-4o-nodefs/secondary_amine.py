"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine has a nitrogen atom bonded to exactly two carbon atoms and one hydrogen atom.

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
        # Check if the atom is nitrogen with hybridization that can fit secondary amines
        if atom.GetAtomicNum() == 7 and atom.GetHybridization() < 4:
            # Get neighbors and verify they include two carbons and one hydrogen
            neighbors = atom.GetNeighbors()
            num_carbon = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 6)
            num_hydrogen = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 1)
            
            # Check the requirements for a secondary amine
            if num_carbon == 2 and num_hydrogen == 1:
                return True, "Contains a secondary amine group"
    
    return False, "Does not contain a secondary amine group"

# You can test the function with a SMILES string to see if it correctly classifies secondary amines.