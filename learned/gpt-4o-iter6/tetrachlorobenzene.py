"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    A tetrachlorobenzene is defined as a molecule containing a benzene ring with
    exactly four chlorine substituents attached at any positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrachlorobenzene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Substructure query for benzene with exactly four chlorine
    benzene_with_chlorine_query = Chem.MolFromSmarts("c1cc(Cl)ccc1Cl")

    # Check for chloro-patterns in the molecule
    if mol.HasSubstructMatch(benzene_with_chlorine_query):
        # Count all chlorine atoms in the structure
        chlorine_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)

        # We need at least four chloro-substituents overall
        if chlorine_count == 4:
            return True, "Benzene ring with exactly four chloro groups found"

    return False, "No benzene ring with exactly four chloro groups found"

# Example usage
smiles = "Clc1cc(Cl)c(Cl)c(Cl)c1"  # Example SMILES for testing
print(is_tetrachlorobenzene(smiles))