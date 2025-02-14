"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is a compound in which a hydroxy group (-OH) is attached to a saturated carbon atom
    which has either three hydrogen atoms attached to it or only one other carbon atom and two hydrogen atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found = False

    for atom in mol.GetAtoms():
        # Look for oxygen atoms with one hydrogen (hydroxyl groups)
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs(includeNeighbors=True) == 1:
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 1:
                continue
            carbon = neighbors[0]
            if carbon.GetAtomicNum() != 6:
                continue
            # Check that carbon is sp3 hybridized (saturated)
            if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Check that carbon is attached to exactly one other carbon (primary carbon)
            num_carbons = sum(1 for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 6)
            if num_carbons != 1:
                continue
            # Check that carbon is attached to exactly two hydrogens
            if carbon.GetTotalNumHs() != 2:
                continue
            found = True
            break

    if not found:
        return False, "No primary alcohol group found"

    return True, "Contains a primary alcohol group attached to a saturated carbon with one other carbon neighbor and two hydrogens"