"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:secondary_amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is characterized by having a nitrogen atom bonded to two hydrocarbyl groups and one hydrogen.

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

    # Iterate over all atoms and find secondary amine nitrogens
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check for nitrogen atoms
            # Count the number of carbon bonds
            carbon_bond_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            non_hydrocarbon_bond_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() != 6)
            
            # Check if there are exactly two carbon bonds
            if carbon_bond_count == 2 and non_hydrocarbon_bond_count <= 1:
                return True, "Contains nitrogen with two hydrocarbyl groups"

    return False, "No nitrogen atom with two hydrocarbyl groups found"