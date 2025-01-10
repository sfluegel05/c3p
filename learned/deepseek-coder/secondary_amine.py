"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:32877 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a secondary amine
    # The pattern looks for a nitrogen atom bonded to exactly two carbon atoms and one hydrogen atom
    secondary_amine_pattern = Chem.MolFromSmarts("[NX3;H1]([#6])([#6])")

    # Check if the molecule contains the secondary amine pattern
    if mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Contains a nitrogen atom bonded to two carbon atoms and one hydrogen atom (secondary amine)"
    else:
        return False, "No secondary amine pattern found"