"""
Classifies: CHEBI:37142 organoiodine compound
"""
"""
Classifies: CHEBI:37143 organoiodine compound
"""
from rdkit import Chem

def is_organoiodine_compound(smiles: str):
    """
    Determines if a molecule is an organoiodine compound based on its SMILES string.
    An organoiodine compound is defined as a compound containing at least one carbon-iodine bond.

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

    # Check for iodine atoms in the molecule
    iodine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53]
    if not iodine_atoms:
        return False, "No iodine atoms found in the molecule"

    # Check if any iodine atom is bonded to at least one carbon atom
    for iodine_atom in iodine_atoms:
        for neighbor in iodine_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon has atomic number 6
                return True, "Contains at least one carbon-iodine bond"

    return False, "No carbon-iodine bonds found"