"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule contains a carboxamidine group based on its SMILES string.
    A carboxamidine has the functional group R-N=C(NR)NR.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a carboxamidine group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a carboxamidine group
    carboxamidine_pattern = Chem.MolFromSmarts("N=C(N)N")
    if carboxamidine_pattern is None:
        return False, "Error in SMARTS pattern"

    # Search for the carboxamidine pattern in the given molecule
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine group"
    else:
        return False, "Does not contain carboxamidine group"