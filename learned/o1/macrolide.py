"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure ring information is computed
    Chem.GetSSSR(mol)

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"

    # Define lactone pattern where both carbonyl carbon and ester oxygen are in a ring
    lactone_smarts = '[C;R](=O)[O;R]'
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)

    # Iterate over rings
    for ring_atoms in atom_rings:
        ring_size = len(ring_atoms)
        if ring_size >= 12:
            # Create a sub-molecule of the ring
            submol = Chem.PathToSubmol(mol, ring_atoms)
            # Search for lactone in the ring sub-molecule
            if submol.HasSubstructMatch(lactone_pattern):
                return True, f"Contains macrocyclic lactone ring of size {ring_size} with lactone functionality"

    return False, "No macrocyclic lactone ring of size 12 or more with lactone functionality found"