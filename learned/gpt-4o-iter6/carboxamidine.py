"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidine has the structure RC(=NR)NR2, where R can be a hydrogen or other substituents.

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

    # Define a generic carboxamidine pattern, allowing for variability in R groups
    carboxamidine_pattern = Chem.MolFromSmarts("C(=N)(N[*])[*]")
    
    # Ensure the molecule has a substructure match with the carboxamidine pattern
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains a carboxamidine structure: RC(=NR)NR2"
    else:
        return False, "Does not contain a carboxamidine structure"