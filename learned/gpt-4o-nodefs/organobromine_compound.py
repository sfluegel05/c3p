"""
Classifies: CHEBI:37141 organobromine compound
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound contains bromine atoms attached to an organic (carbon-containing) structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for bromine atoms
    has_bromine = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 35:  # Atomic number for bromine
            # Check if the bromine atom is connected to at least one carbon atom
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Atomic number for carbon
                    has_bromine = True
                    break
            if has_bromine:
                break

    if has_bromine:
        return True, "Contains bromine atoms attached to a carbon-based structure"
    else:
        return False, "No bromine atoms attached to an organic structure found"