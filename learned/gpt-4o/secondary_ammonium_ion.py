"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is characterized by a nitrogen atom with a +1 charge attached
    to two distinct carbon chains, typically formed by protonation of a secondary amine (R2NH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for secondary ammonium ion
    # Specifically targeting secondary amines with a formal +1 charge and precisely two carbon atom bonds
    secondary_ammonium_pattern = Chem.MolFromSmarts("[N;+1]([C;H3,H2])([C;H3,H2])")  # Nitrogen with +1, attached to two non-aromatic carbons

    if mol.HasSubstructMatch(secondary_ammonium_pattern):
        return True, "Contains protonated secondary ammonium ion pattern"

    return False, "Does not contain a secondary ammonium ion pattern"