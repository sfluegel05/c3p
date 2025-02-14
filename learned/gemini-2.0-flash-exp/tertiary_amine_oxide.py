"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine oxide based on its SMILES string.
    A tertiary amine oxide is an N-oxide where there are three organic groups bonded to the nitrogen atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine oxide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for N-oxide using SMARTS
    n_oxide_pattern = Chem.MolFromSmarts("[N+][O-]")
    if not mol.HasSubstructMatch(n_oxide_pattern):
        return False, "Not an N-oxide"

    # Find the nitrogen atom and oxygen atom from the N-oxide bond.
    match = mol.GetSubstructMatch(n_oxide_pattern)
    nitrogen_atom = mol.GetAtomWithIdx(match[0])
    oxygen_atom = mol.GetAtomWithIdx(match[1])

    # Check that nitrogen has 4 bonds (3 with carbon, 1 with oxygen)
    if nitrogen_atom.GetTotalValence() != 4:
        return False, f"Nitrogen has {nitrogen_atom.GetTotalValence()} connections, not 4"

    # Check the connectivity.
    neighbors = nitrogen_atom.GetNeighbors()

    # Ensure the neighbors bonded to nitrogen are carbon, and not the oxygen which has already been verified.
    carbon_neighbors = []
    for neighbor in neighbors:
        if neighbor.GetIdx() != oxygen_atom.GetIdx(): # Ensure we are not evaluating the oxygen again.
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors.append(neighbor)
            else:
                return False, "Nitrogen not bonded to three carbon atoms."
    if len(carbon_neighbors) != 3:
        return False, "Nitrogen not bonded to three carbon atoms."

    # Check that the carbon neighbours are connected to a complex organic group, and not just hydrogen atoms.
    for carbon_neighbor in carbon_neighbors:
        neighbor_atoms = carbon_neighbor.GetNeighbors()
        is_organic = False
        for neighbor_atom in neighbor_atoms:
            if neighbor_atom.GetIdx() != nitrogen_atom.GetIdx(): # Avoid checking back towards the central N atom
                if neighbor_atom.GetAtomicNum() != 1:
                    is_organic = True
        if not is_organic:
             return False, "Nitrogen is not bonded to three organic groups"

    return True, "Tertiary amine oxide detected"