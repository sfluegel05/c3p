"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:82795 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is a benzene ring with exactly four chlorine atoms attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzene ring pattern
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(benzene_pattern):
        return False, "No benzene ring found"

    # Count the number of chlorine atoms attached to the benzene ring
    chlorine_count = 0
    benzene_atoms = mol.GetSubstructMatch(benzene_pattern)
    for atom_idx in benzene_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 17:  # Chlorine atomic number
                chlorine_count += 1

    if chlorine_count != 4:
        return False, f"Found {chlorine_count} chlorine atoms, need exactly 4"

    return True, "Contains a benzene ring with exactly four chlorine atoms attached"