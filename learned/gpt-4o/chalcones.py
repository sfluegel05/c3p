"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    Chalcones are defined as 1,3-diphenylpropenone or benzylideneacetophenone derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a chalcone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more focused pattern for the chalcone core ArCH=CHC(=O)Ar
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1/C=C/C(=O)c2ccccc2")
    
    # Attempt to match the chalcone pattern in the molecule
    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Molecule matches chalcone core structure ArCH=CHC(=O)Ar"
        
    return False, "Does not match chalcone core structure"