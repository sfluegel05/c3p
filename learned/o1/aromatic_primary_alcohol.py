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

    found = False  # Flag to indicate if aromatic primary alcohol is found

    # Iterate over all atoms to find hydroxyl groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            # Check if oxygen is part of a hydroxyl group (bonded to hydrogen)
            has_hydrogen = any(nb.GetAtomicNum() == 1 for nb in atom.GetNeighbors())
            if not has_hydrogen:
                continue  # Not a hydroxyl group
            # Get the carbon atom connected to this oxygen
            neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 6]
            if len(neighbors) != 1:
                continue  # Oxygen not connected to a carbon (unlikely but check)
            carbon = neighbors[0]
            # Check if carbon is sp3 hybridized (indicative of single bonds)
            if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue  # Carbon is not sp3 hybridized
            # Check if carbon has exactly two hydrogens (primary alcohol)
            h_count = sum(1 for nb in carbon.GetNeighbors() if nb.GetAtomicNum() == 1)
            if h_count != 2:
                continue  # Not a primary alcohol
            # Check if carbon is bonded to an aromatic ring
            aromatic_bonded = any(
                nb.GetIsAromatic() and nb.GetAtomicNum() in [6,7,8,15,16]
                for nb in carbon.GetNeighbors()
                if nb.GetIdx() != atom.GetIdx() and nb.GetAtomicNum() != 1
            )
            if aromatic_bonded:
                found = True
                break  # No need to check further

    if found:
        return True, "Molecule is an aromatic primary alcohol"
    else:
        return False, "No aromatic primary alcohol group found"