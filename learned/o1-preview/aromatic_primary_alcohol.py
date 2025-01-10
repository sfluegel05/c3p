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
    is attached to a carbon which is itself bonded to an aromatic carbon atom in a ring.

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

    # Identify primary alcohol carbons (sp3 carbon with two hydrogens bonded to OH group)
    primary_alcohol_pattern = Chem.MolFromSmarts('[C;H2;X4][OX2H]')
    matches = mol.GetSubstructMatches(primary_alcohol_pattern)

    if not matches:
        return False, "No primary alcohol groups found"

    # Iterate over all primary alcohol carbons
    for match in matches:
        carbon_idx = match[0]  # Index of carbon atom bearing OH
        oxygen_idx = match[1]  # Index of oxygen atom

        carbon_atom = mol.GetAtomWithIdx(carbon_idx)

        # Check neighbors excluding the oxygen atom
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetIdx() == oxygen_idx:
                continue  # Skip the oxygen atom

            # Check if neighbor is an aromatic carbon atom in a ring
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic() and neighbor.IsInRing():
                return True, "Molecule contains an aromatic primary alcohol group"

    # No matching groups found
    return False, "No aromatic primary alcohol groups connected to aromatic carbon atom found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:33853',
        'name': 'aromatic primary alcohol',
        'definition': 'Any primary alcohol in which the alcoholic hydroxy group is attached to a carbon which is itself bonded to an aromatic ring.',
        'parents': ['CHEBI:33824', 'CHEBI:33854']
    }
}