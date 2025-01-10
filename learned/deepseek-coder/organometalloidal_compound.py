"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: CHEBI:50860 organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound is defined as a compound having bonds between one or more
    metalloid atoms and one or more carbon atoms of an organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define metalloid atoms (As, Si, B, Ge, Sb, Te, Po)
    metalloid_atomic_numbers = {33, 14, 5, 32, 51, 52, 84}

    # Check if the molecule contains any metalloid atoms
    metalloid_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() in metalloid_atomic_numbers]
    if not metalloid_atoms:
        return False, "No metalloid atoms found"

    # Check for at least one carbon-metalloid bond where the carbon is part of an organyl group
    for atom in metalloid_atoms:
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                # Ensure the carbon is part of an organyl group (i.e., bonded to at least one hydrogen or another carbon)
                if any(nbr.GetAtomicNum() == 1 or nbr.GetAtomicNum() == 6 for nbr in neighbor.GetNeighbors()):
                    return True, "Contains at least one carbon-metalloid bond in an organyl group"

    return False, "No carbon-metalloid bonds found in an organyl group"