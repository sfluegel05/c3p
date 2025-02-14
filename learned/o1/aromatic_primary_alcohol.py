"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:33854 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol is a primary alcohol where the alcoholic hydroxy group is attached
    to a carbon which is directly bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to accurately count hydrogen atoms
    mol = Chem.AddHs(mol)

    found = False  # Flag to indicate if the pattern is found

    for atom in mol.GetAtoms():
        # Step 1: Identify hydroxyl groups (-OH)
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue

            # Check if oxygen is connected to one hydrogen and one carbon
            attached_hydrogen = None
            alcoholic_carbon = None
            for neighbor in neighbors:
                if neighbor.GetAtomicNum() == 1:
                    attached_hydrogen = neighbor
                elif neighbor.GetAtomicNum() == 6:
                    alcoholic_carbon = neighbor

            if attached_hydrogen is None or alcoholic_carbon is None:
                continue

            # Step 2: Verify that the alcoholic carbon is primary
            carbon_neighbors = [n for n in alcoholic_carbon.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() != atom.GetIdx()]
            if len(carbon_neighbors) != 1:
                continue  # Not a primary carbon

            # Step 3: Check if the alcoholic carbon is connected to an aromatic ring
            is_connected_to_aromatic = False
            for neighbor in carbon_neighbors:
                if neighbor.GetIsAromatic():
                    is_connected_to_aromatic = True
                    break

            if is_connected_to_aromatic:
                found = True
                break  # Pattern found

    if found:
        return True, "Molecule is an aromatic primary alcohol"
    else:
        return False, "Does not meet criteria for an aromatic primary alcohol"