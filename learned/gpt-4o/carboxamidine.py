"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine is characterized by the presence of a -C(=NR)NR2 group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general carboxamidine SMARTS pattern
    # -C(=NR)NR2, this considers any carbon (can be aromatic), with a double bonded N and another N
    carboxamidine_pattern = Chem.MolFromSmarts("[#6]=N[#7]")

    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine group (-C(=NR)NR2)"

    return False, "Does not contain a carboxamidine group"

# Example usage:
# test_smiles = "NC(=N)c1ccccc1"  # benzamidine
# print(is_carboxamidine(test_smiles))