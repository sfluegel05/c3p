"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine is a compound formally derived from ammonia by replacing one, two or three hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for nitrogen atoms in the molecule
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found in the molecule"

    # Check that nitrogen atoms have hydrocarbyl groups
    for atom in nitrogen_atoms:
        # Get neighbors of the nitrogen atom
        neighbors = atom.GetNeighbors()
        # Count carbon neighbors
        carbon_neighbors = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 6)
        # Check for amine type
        if carbon_neighbors >= 1:
            return True, f"Contains nitrogen with {carbon_neighbors} hydrocarbyl groups suggesting it is an amine"
    
    return False, "Nitrogen is not bonded to hydrocarbyl groups as required for amines"