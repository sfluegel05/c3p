"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone or its derivatives based on its SMILES string.
    A chalcone is characterized by a 1,3-diphenylpropenone structure (ArCH=CH(=O)Ar)
    that may include various aromatic derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Primary chalcone and loosely substituted form: Ar-C=CH-C(=O)-Ar
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1C=CC(=O)C2=CC=CC=C2")
    substituted_chalcone_pattern = Chem.MolFromSmarts("C:c(:c):C=C:C(=O):c(:c):c")

    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Contains the primary 1,3-diphenylpropenone structure"
    elif mol.HasSubstructMatch(substituted_chalcone_pattern):
        return True, "Contains a variant structure with substituted aromatic rings on chalcone backbone"

    return False, "Does not contain a chalcone structure"