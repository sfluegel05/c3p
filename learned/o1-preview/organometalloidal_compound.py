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
    An organometalloidal compound is defined as 'A compound having bonds between one or more metalloid atoms
    and one or more carbon atoms of an organyl group.'

    Metalloid elements considered here are arsenic (As) and antimony (Sb).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    # List of metalloid atomic numbers to consider
    metalloid_atomic_nums = [33, 51]  # As, Sb

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over metalloid atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metalloid_atomic_nums:
            metalloid_atom = atom
            has_carbon_neighbor = False
            has_oxygen_neighbor = False
            # Check neighbors of metalloid atom
            for neighbor in metalloid_atom.GetNeighbors():
                atomic_num = neighbor.GetAtomicNum()
                if atomic_num == 6:  # Carbon atom
                    has_carbon_neighbor = True
                elif atomic_num == 8:  # Oxygen atom
                    has_oxygen_neighbor = True
            if has_carbon_neighbor and has_oxygen_neighbor:
                return True, f"Metalloid atom ({metalloid_atom.GetSymbol()}) bonded to carbon and oxygen atoms"
    # If no metalloid-carbon bonds found
    return False, "No metalloid atom bonded to both carbon and oxygen atoms of organyl groups found"