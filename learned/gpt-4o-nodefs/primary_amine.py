"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine contains a nitrogen atom that is bonded to exactly one carbon 
    and two hydrogen atoms (RNH2).

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

    # Iterate through atoms to find potential primary amine groups
    for atom in mol.GetAtoms():
        # Look for nitrogen atoms
        if atom.GetAtomicNum() == 7:
            # Primary amine Nitrogen should have exactly three neighbors
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                # Count carbon and hydrogen neighbors
                carbon_count = sum(1 for nbr in neighbors if nbr.GetAtomicNum() == 6)
                hydrogen_count = sum(1 for nbr in neighbors if nbr.GetAtomicNum() == 1)
                # A primary amine should have one carbon and two hydrogens bonded to nitrogen
                if carbon_count == 1 and hydrogen_count == 2:
                    return True, "Contains a primary amine group (one carbon and two hydrogens bonded to nitrogen)"
    
    return False, "No primary amine group found"