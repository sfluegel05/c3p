"""
Classifies: CHEBI:37142 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound is characterized by one or more iodine atoms bonded
    directly to a carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organoiodine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for iodine atoms within the molecule
    has_iodine = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 53:  # Atomic number of iodine is 53
            # Check if the iodine is bonded directly to a carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Atomic number of carbon is 6
                    has_iodine = True
                    break
            if has_iodine:
                break

    if has_iodine:
        return True, "Iodine is bonded directly to carbon, indicating an organoiodine compound."

    return False, "No iodine-carbon bonds found; not an organoiodine compound."