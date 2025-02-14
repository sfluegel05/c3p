"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: methyl sulfide
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide where the sulfur atom is connected to two carbon atoms,
    and at least one of those carbons is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize flag
    is_methyl_sulfide = False

    # Identify sulfur atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            # Check if sulfur is attached to exactly two carbon atoms
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 2:
                continue
            carbon_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 2:
                continue
            # Check if at least one neighboring carbon is a methyl group
            for carbon in carbon_neighbors:
                if carbon.GetDegree() == 1:
                    # Carbon attached only to sulfur (methyl group)
                    return True, "Sulfur atom attached to a methyl group and another carbon"
    return False, "No sulfur atom attached to a methyl group and another carbon"