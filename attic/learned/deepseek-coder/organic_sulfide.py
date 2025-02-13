"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: CHEBI:48376 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide is defined as a compound with the structure R-S-R' (R, R' =/= H).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sulfur atoms bonded to two carbon atoms (R-S-R')
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]  # Sulfur has atomic number 16
    for sulfur in sulfur_atoms:
        neighbors = sulfur.GetNeighbors()
        if len(neighbors) == 2:  # Sulfur must be bonded to exactly two atoms
            if all(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):  # Both neighbors must be carbon atoms
                return True, "Contains R-S-R' pattern (R, R' =/= H)"

    return False, "No R-S-R' pattern found (R, R' =/= H)"