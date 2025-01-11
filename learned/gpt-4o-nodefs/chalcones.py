"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone has two aromatic rings connected by a three-carbon α,β-unsaturated carbonyl group (C=O-C=C).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for chalcones
    chalcone_pattern = Chem.MolFromSmarts("c1ccccc1C(=O)C=Cc2ccccc2")

    # Check for chalcone pattern in the molecule
    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Contains chalcone structure: two aromatic rings connected by α,β-unsaturated carbonyl"

    return False, "Does not contain chalcone structure: missing characteristic α,β-unsaturated carbonyl linkage"