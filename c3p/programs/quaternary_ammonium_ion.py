"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion based on its SMILES string.
    A quaternary ammonium ion has a central nitrogen atom with four organic substituents and a positive charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a nitrogen with four carbon attachments and a positive charge
    # Each substituent must be a carbon atom (including aromatic or aliphatic)
    quaternary_n_pattern = Chem.MolFromSmarts("[N+]([#6])([#6])([#6])[#6]")
    
    # Check if the molecule contains the pattern
    if mol.HasSubstructMatch(quaternary_n_pattern):
        return True, "Contains a nitrogen atom with four organic substituents and a positive charge"
    
    return False, "No quaternary ammonium group found"