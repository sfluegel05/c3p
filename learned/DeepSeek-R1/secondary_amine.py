"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:48377 secondary amine
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

    # Iterate through all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Check for exactly three bonds (two carbons + one hydrogen)
            if atom.GetDegree() == 3:
                carbon_neighbors = 0
                # Count carbon neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Carbon atom
                        carbon_neighbors += 1
                # Secondary amine requires exactly two carbon neighbors
                if carbon_neighbors == 2:
                    return True, "Contains a nitrogen atom with exactly two carbon neighbors and three bonds"
    
    return False, "No secondary amine group (N with two carbon substituents) found"