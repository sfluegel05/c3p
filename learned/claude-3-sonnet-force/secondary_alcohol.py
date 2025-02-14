"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: CHEBI:35569 secondary alcohol

A secondary alcohol is a compound in which a hydroxy group, -OH, is attached to a saturated carbon atom which has two or more other carbon atoms attached to it.
"""

from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary alcohol group, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find atoms with -OH group
    oh_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and sum(mol.GetAtomWithIdx(n).GetTotalNumHs() for n in atom.GetNeighbors()) == 1]

    # Check if any of those atoms are bound to a saturated carbon with at least 2 other carbons
    for oh_atom in oh_atoms:
        carbon = mol.GetAtomWithIdx(list(oh_atom.GetNeighbors())[0])
        if carbon.GetAtomicNum() == 6 and carbon.GetTotalNumHs() + len(carbon.GetNeighbors()) >= 4 and sum(1 for n in carbon.GetNeighbors() if mol.GetAtomWithIdx(n).GetAtomicNum() == 6) >= 2:
            return True, "Contains a hydroxy group attached to a saturated carbon with at least two other carbon atoms"

    return False, "No secondary alcohol group found"