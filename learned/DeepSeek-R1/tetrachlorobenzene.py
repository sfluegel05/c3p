"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is a chlorobenzene with exactly four chlorine substituents on a benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all benzene rings (aromatic six-membered carbon rings)
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)
    if not benzene_matches:
        return False, "No benzene ring found"
    
    # Check each benzene ring for exactly four chlorine substituents
    for ring in benzene_matches:
        cl_count = 0
        # Iterate through each atom in the benzene ring
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check all neighbors of the ring atom
            for neighbor in atom.GetNeighbors():
                # Count chlorine atoms not part of the ring
                if neighbor.GetAtomicNum() == 17 and neighbor.GetIdx() not in ring:
                    cl_count += 1
        if cl_count == 4:
            return True, "Contains a benzene ring with four chlorine substituents"
    
    return False, "No benzene ring found with exactly four chlorine substituents"