"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is characterized by a nitrogen atom bonded to exactly 
    three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for tertiary amine pattern
    is_tertiary_amine = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check if the atom is nitrogen
            bonded_carbon_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6)
            if bonded_carbon_count == 3:  # Nitrogen bonded to exactly three carbon atoms
                is_tertiary_amine = True
                break

    if is_tertiary_amine:
        return True, "Molecule contains a tertiary amine (N bonded to 3 carbon atoms)"
    else:
        return False, "No tertiary amine found (N not bonded to 3 carbon atoms)"