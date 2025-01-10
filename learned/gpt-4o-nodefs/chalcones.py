"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone has two aromatic rings connected by a three-carbon α,β-unsaturated carbonyl group.

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

    # Define a more flexible SMARTS pattern for chalcones
    # Look for two aromatic rings connected by an α,β-unsaturated carbonyl group
    aromatic_ring_pattern = "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1"
    enone_pattern = "[#6](=O)C=C"
    chalcone_pattern = aromatic_ring_pattern + enone_pattern + aromatic_ring_pattern
    
    # Convert SMARTS strings to pattern objects
    generic_chalcone_pattern = Chem.MolFromSmarts(chalcone_pattern)
    
    # Check for chalcone pattern in the molecule with flexibility for substitutions
    if mol.HasSubstructMatch(generic_chalcone_pattern):
        return True, "Contains chalcone structure: two aromatic rings connected by α,β-unsaturated carbonyl"

    return False, "Does not contain chalcone structure: missing characteristic α,β-unsaturated carbonyl linkage"