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
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for primary alcohol carbon (CH2 bonded to OH)
    primary_alcohol_pattern = Chem.MolFromSmarts("[C;H2;X4][OH]")
    matches = mol.GetSubstructMatches(primary_alcohol_pattern)

    if not matches:
        return False, "No primary alcohol groups found"

    # Iterate over all primary alcohol groups found
    for match in matches:
        carbon_idx = match[0]  # Index of the carbon atom
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Check if the carbon atom is connected to an aromatic ring
        is_connected_to_aromatic = False
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                is_connected_to_aromatic = True
                break  # No need to check further if one aromatic neighbor is found
        
        if is_connected_to_aromatic:
            return True, "Molecule contains an aromatic primary alcohol group"

    # If no primary alcohol carbons are connected to an aromatic ring
    return False, "No aromatic primary alcohol groups connected to aromatic ring found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33853',
        'name': 'aromatic primary alcohol',
        'definition': 'Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon which is itself bonded to an aromatic ring.',
        'parents': ['CHEBI:33824', 'CHEBI:33854']
    }
}