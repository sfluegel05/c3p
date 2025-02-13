"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine has a nitrogen atom connected to two carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Atomic number of nitrogen
            if atom.GetTotalDegree() == 3:  # Total degree includes implicit hydrogens
                carbon_neighbor_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
                if carbon_neighbor_count == 2:
                    return True, "Nitrogen is bonded to two carbons and has a total degree of 3, indicative of a secondary amine"
    
    return False, "No nitrogen atom bonded to exactly two carbons found"