"""
Classifies: CHEBI:139588 alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:35551 alpha-hydroxy ketone
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is an alpha-hydroxy ketone based on its SMILES string.
    An alpha-hydroxy ketone is a ketone containing a hydroxy group on the alpha-carbon relative to the C=O group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carbonyl atoms
    carbonyl_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 0]

    # Check each carbonyl for alpha-hydroxy group
    for carbonyl in carbonyl_atoms:
        alpha_carbon = None
        for neighbor in carbonyl.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbon = neighbor
                break
        if alpha_carbon is None:
            continue

        # Check if alpha carbon has a hydroxy group
        has_hydroxy = False
        for neighbor in alpha_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
                has_hydroxy = True
                break

        if has_hydroxy:
            return True, "Contains a carbonyl group with a hydroxy group on the alpha-carbon"

    return False, "No alpha-hydroxy ketone group found"