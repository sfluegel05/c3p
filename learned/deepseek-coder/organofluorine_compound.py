"""
Classifies: CHEBI:37143 organofluorine compound
"""
from rdkit import Chem

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.
    An organofluorine compound is defined as a compound containing at least one carbon-fluorine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for any fluorine atoms in the molecule
    fluorine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9]
    if not fluorine_atoms:
        return False, "No fluorine atoms found"

    # Check if any fluorine atom is bonded to a carbon atom
    for atom in fluorine_atoms:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                return True, "Contains at least one carbon-fluorine bond"

    return False, "No carbon-fluorine bonds found"