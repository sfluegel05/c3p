"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: CHEBI:39277 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of bromine atoms
    bromine_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 35]  # Atomic number of bromine is 35
    if not bromine_atoms:
        return False, "No bromine atoms found"

    # Check if any bromine atom is directly bonded to a carbon atom
    for bromine_atom in bromine_atoms:
        for neighbor in bromine_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Atomic number of carbon is 6
                return True, "Contains at least one carbon-bromine bond"

    return False, "No carbon-bromine bonds found"