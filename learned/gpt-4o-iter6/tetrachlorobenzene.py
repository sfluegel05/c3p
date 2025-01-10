"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a molecule containing an aromatic ring with exactly four chloro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure query for any aromatic ring
    aromatic_ring_query = Chem.MolFromSmarts("a")
    chloro_count = 0

    # Iterate over the flags set on each atom to determine if they are
    # part of an aromatic system and if so, count the chlorines within
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            # Collect atoms in the same aromatic ring
            ring_atoms = [a.GetIdx() for a in atom.GetNeighbors() if a.GetIsAromatic()]

            # Count chlorine atoms on the aromatic system
            for ring_atom in ring_atoms:
                neighbor = mol.GetAtomWithIdx(ring_atom)
                if neighbor.GetAtomicNum() == 17:  # Atomic number for Cl
                    chloro_count += 1

            # If there are exactly four chlorines on this aromatic ring, classify
            if chloro_count == 4:
                return True, "Aromatic ring with exactly four chloro groups found"

    # If no suitable aromatic ring is found, return False
    return False, "No aromatic ring with exactly four chloro groups found"

# Example usage
smiles = "Clc1ccc(Cl)c(Cl)c1Cl"  # Example SMILES for testing
print(is_tetrachlorobenzene(smiles))