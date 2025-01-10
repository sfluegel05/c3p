"""
Classifies: CHEBI:29256 thiol
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group, -SH, is attached to a carbon atom of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule
    mol_with_H = Chem.AddHs(mol)

    # Search for thiol groups (-SH attached to carbon)
    thiol_found = False
    for atom in mol_with_H.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Get neighbors of sulfur atom
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue  # Sulfur must be connected to exactly two atoms
            hydrogen_neighbor = False
            carbon_neighbor = False
            for neighbor in neighbors:
                atomic_num = neighbor.GetAtomicNum()
                if atomic_num == 1:
                    hydrogen_neighbor = True
                elif atomic_num == 6:
                    carbon_neighbor = True
            if hydrogen_neighbor and carbon_neighbor:
                thiol_found = True
                break

    if thiol_found:
        return True, "Contains thiol group (-SH) attached to carbon atom"
    else:
        return False, "No thiol group (-SH) attached to carbon atom found"