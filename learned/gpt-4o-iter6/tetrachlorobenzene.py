"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a benzene ring carrying four chloro groups at unspecified positions.

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

    # Find benzene rings in the molecule
    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
    benzene_rings = mol.GetSubstructMatches(benzene_pattern)
    
    if not benzene_rings:
        return False, "No benzene rings found"

    # Check for chlorine atoms attached to each benzene ring
    for ring in benzene_rings:
        chlorine_count = sum(1 for atom_idx in ring if mol.GetAtomWithIdx(atom_idx).GetSymbol() == 'Cl')
        
        # Count chlorines directly connected to the benzene ring
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'Cl' and neighbor.GetIdx() not in ring:
                    chlorine_count += 1
        
        if chlorine_count == 4:
            return True, "Benzene ring with exactly four chloro groups found"
        
    return False, "No benzene ring with exactly four chloro groups found"