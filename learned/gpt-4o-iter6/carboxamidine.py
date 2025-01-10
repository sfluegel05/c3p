"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidine has the structure RC(=NR)NR2, denoting the -C(=NH)NH2 group.

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

    # Define the carboxamidine pattern
    carboxamidine_pattern = Chem.MolFromSmarts("N=C(N)N")
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains the carboxamidine group: N=C(N)N"
    else:
        return False, "Does not contain the carboxamidine group"