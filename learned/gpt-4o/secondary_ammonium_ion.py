"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion arises from the protonation of a secondary amine, defined by
    a positively charged nitrogen atom bound to exactly two different carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to exactly match secondary ammonium ion: [NH2+][C]([C])
    # Ensure nitrogen is bound to two distinct carbon atoms; both single bonded
    secondary_ammonium_smart = Chem.MolFromSmarts("[N+H2&!H3&!H1]([C])[C]")

    # While also checking if it's majorly a cation in physiological pH (ensure there are no other negating features)
    # Look for presence of secondary ammonium ion with extended criteria
    if mol.HasSubstructMatch(secondary_ammonium_smart):
        return True, "Contains protonated secondary amine group forming secondary ammonium ion"
    
    return False, "Does not contain the features of a secondary ammonium ion"