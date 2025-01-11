"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is defined as a nitrogen atom connected to three hydrocarbyl groups, which in SMILES terms means a nitrogen bonded to three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string to obtain a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify nitrogen atoms in the molecule
    is_tertiary_amine_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check for nitrogen
            # Count the number of carbon bonds
            carbon_count = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            # Ensure the nitrogen is bonded to exactly three carbon atoms
            if carbon_count == 3:
                is_tertiary_amine_found = True
                break
    
    if is_tertiary_amine_found:
        return True, "Nitrogen bonded to three carbon atoms found indicating a tertiary amine"
    else:
        return False, "No nitrogen bonded to three carbon atoms found"