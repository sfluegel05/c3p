"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:35671 chalcones
"""

from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is an aromatic ketone that consists of two aromatic rings
    joined by a three-carbon α,β-unsaturated carbonyl system (Ar-CH=CH-C(=O)-Ar),
    and its derivatives formed by substitution.

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

    # Define the chalcone core SMARTS pattern
    # Pattern represents:
    # Aromatic ring connected to C=C-C=O connected to aromatic ring
    chalcone_pattern = Chem.MolFromSmarts("[a][C]=[C]-C(=O)[a]")
    if chalcone_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for the chalcone substructure
    if mol.HasSubstructMatch(chalcone_pattern):
        return True, "Contains chalcone core structure"
    else:
        return False, "Does not contain chalcone core structure"