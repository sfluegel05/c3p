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
    An organometalloidal compound is defined as 'A compound having bonds between one or more metalloid atoms and one or more carbon atoms of an organyl group.'

    Metalloid elements include boron (B), silicon (Si), germanium (Ge), arsenic (As), antimony (Sb), and tellurium (Te).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    # List of metalloid atomic numbers
    metalloid_atomic_nums = [5, 14, 32, 33, 51, 52]  # B, Si, Ge, As, Sb, Te

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flag to check for metalloid-carbon bond
    has_metalloid_carbon_bond = False

    # Iterate over atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() in metalloid_atomic_nums:
            # This is a metalloid atom
            metalloid_atom = atom
            # Get neighbors
            neighbors = atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 6:
                    # Metalloid atom is bonded to carbon
                    return True, f"Metalloid atom ({atom.GetSymbol()}) bonded to carbon atom"
    # If no metalloid-carbon bonds found
    return False, "No metalloid-carbon bonds found"