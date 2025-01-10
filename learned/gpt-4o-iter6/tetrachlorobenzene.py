"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a molecule containing one or more benzene rings,
    where one of the rings has exactly four chlorine substituents attached at any positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Substructure query for aromatic ring with four chlorine substitutions
    # This pattern allows for arbitrary positioning on benzene rings 
    benzene_query = Chem.MolFromSmarts("c1ccccc1")
    all_benzene_rings = mol.GetSubstructMatches(benzene_query)
    
    # Iterate over each benzene ring detected
    for ring in all_benzene_rings:
        chlorine_attached_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check neighbors of aromatic carbon atoms in the ring
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 17:  # Atomic number of Chlorine is 17
                    chlorine_attached_count += 1

        # Check for exactly four chlorine atoms on this aromatic ring
        if chlorine_attached_count == 4:
            return True, "Benzene ring with exactly four chloro groups found"
    
    return False, "No benzene ring with exactly four chloro groups found"

# Example usage
smiles = "Clc1cc(Cl)c(Cl)c(Cl)c1"  # Example SMILES for testing
print(is_tetrachlorobenzene(smiles))