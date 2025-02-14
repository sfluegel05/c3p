"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to exactly two carbon atoms and one hydrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            if atom.GetFormalCharge() != 0:
                continue  # Skip if nitrogen has a formal charge
            total_degree = atom.GetTotalDegree()  # Includes hydrogens
            if total_degree != 3:
                continue  # Skip if nitrogen is not connected to exactly 3 atoms
            h_count = atom.GetTotalNumHs()
            if h_count != 1:
                continue  # Skip if nitrogen doesn't have exactly one hydrogen
            # Count the number of carbon neighbors
            c_count = 0
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Carbon atom
                    c_count += 1
            if c_count == 2:
                return True, "Contains secondary amine group"
    return False, "Does not contain secondary amine group"