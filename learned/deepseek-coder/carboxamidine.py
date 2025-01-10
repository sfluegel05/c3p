"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:37671 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine has the structure RC(=NR)NR2.

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

    # Define the carboxamidine pattern: RC(=NR)NR2
    carboxamidine_pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3H2,NX3H1,NX3H0]")
    if not mol.HasSubstructMatch(carboxamidine_pattern):
        return False, "No carboxamidine structure (RC(=NR)NR2) found"

    return True, "Contains the carboxamidine structure (RC(=NR)NR2)"