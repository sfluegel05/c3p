"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidine has the structure RC(=NR)NR2, denoting the -C(=NR)NR2 group,
    where R can be hydrogen or other substituents, such as alkyl or aryl groups.

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

    # Define a more general carboxamidine pattern
    carboxamidine_pattern = Chem.MolFromSmarts("C(=N)(N)N")
    
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains a carboxamidine-like structure: C(=N)(N)N"
    else:
        return False, "Does not contain a carboxamidine-like structure"