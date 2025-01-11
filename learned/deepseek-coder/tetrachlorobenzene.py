"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: CHEBI:82795 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule contains a tetrachlorobenzene moiety based on its SMILES string.
    A tetrachlorobenzene moiety is a benzene ring with exactly four chlorine atoms attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a tetrachlorobenzene moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all benzene rings in the molecule
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_matches:
        return False, "No benzene ring found"

    # Check each benzene ring for exactly four chlorine substituents
    for match in benzene_matches:
        # Get atoms in this benzene ring
        ring_atoms = set(match)
        
        # Count chlorine atoms directly attached to this benzene ring
        chlorine_count = 0
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17 and neighbor.GetIdx() not in ring_atoms:
                    chlorine_count += 1
        
        if chlorine_count == 4:
            return True, "Contains a benzene ring with exactly four chlorine atoms attached"

    return False, "No benzene ring with exactly four chlorine atoms found"