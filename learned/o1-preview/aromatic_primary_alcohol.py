"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:33853 aromatic primary alcohol
Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon which is itself bonded to an aromatic ring.
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.

    An aromatic primary alcohol is defined as any primary alcohol in which the alcoholic hydroxy group 
    is attached to a carbon which is itself bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an aromatic primary alcohol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern to find all carbons attached to an OH group
    alcohol_pattern = Chem.MolFromSmarts("[C][OH]")
    matches = mol.GetSubstructMatches(alcohol_pattern)

    if not matches:
        return False, "No alcohol groups found"

    # Iterate over all alcohol carbons
    for match in matches:
        carbon_idx = match[0]  # Index of carbon atom bearing OH
        oxygen_idx = match[1]  # Index of oxygen atom

        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # Check if the carbon atom is primary (connected to only one other carbon)
        num_carbon_neighbors = 0
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbon_idx:
                num_carbon_neighbors += 1
        if num_carbon_neighbors != 1:
            continue  # Not a primary alcohol

        # Check if the carbon atom is directly connected to an aromatic ring
        is_connected_to_aromatic = False
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                is_connected_to_aromatic = True
                break
        if not is_connected_to_aromatic:
            continue  # Not connected to aromatic ring

        # Passed all checks
        return True, "Molecule contains an aromatic primary alcohol group"

    # No matching groups found
    return False, "No aromatic primary alcohol groups connected to aromatic ring found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33853',
        'name': 'aromatic primary alcohol',
        'definition': 'Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon which is itself bonded to an aromatic ring.',
        'parents': ['CHEBI:33824', 'CHEBI:33854']
    }
}