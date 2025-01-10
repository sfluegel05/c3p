"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between metalloid atoms (e.g., arsenic) 
    and carbon atoms of an organyl group.

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
    
    # Define metalloid atomic numbers, focusing on arsenic for current examples
    metalloid_atomic_number = 33  # Atomic number for arsenic

    # Check for the presence of a metalloid atom (arsenic) directly bonded to carbon
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == metalloid_atomic_number:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:  # Check if neighbor is a carbon
                    return True, "Contains metalloid (arsenic) directly bonded to carbon: an organometalloidal compound"
    
    return False, "No suitable metalloid-carbon bonds identified"