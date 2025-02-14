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
    to a carbon which is itself bonded directly to an aromatic ring.

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

    # Define SMARTS pattern for primary alcohol connected to aromatic ring
    # [CX3H2]-[OH] represents a primary alcohol carbon
    # [CX3H2]-[OH]-[*]~[*]: The alcoholic carbon is connected to any atom (*) that is part of an aromatic ring (~[*]: any bond to an aromatic atom)
    pattern = Chem.MolFromSmarts('[CX3H2][OH]')

    matches = mol.GetSubstructMatches(pattern)
    if not matches:
        return False, "No primary alcohol group found"

    for match in matches:
        alcohol_carbon_idx = match[0]  # Index of the alcoholic carbon
        alcohol_carbon = mol.GetAtomWithIdx(alcohol_carbon_idx)

        # Check if the alcoholic carbon is connected to an aromatic ring
        aromatic_neighbor = False
        for neighbor in alcohol_carbon.GetNeighbors():
            if neighbor.GetIdx() != match[1]:  # Exclude the oxygen atom in the hydroxyl group
                if neighbor.GetIsAromatic():
                    # Further check if neighbor is part of an aromatic ring
                    rings = neighbor.GetOwningMol().GetRingInfo()
                    if rings.IsAtomInRingOfSize(neighbor.GetIdx(), 6) or rings.IsAtomInRingOfSize(neighbor.GetIdx(), 5):
                        aromatic_neighbor = True
                        break

        if aromatic_neighbor:
            return True, "Molecule is an aromatic primary alcohol"

    return False, "No aromatic primary alcohol group found"