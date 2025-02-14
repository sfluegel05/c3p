"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as an organic compound with two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for amino group (nitrogen with 0, 1 or 2 hydrogens) and 3 bonds
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0]")
    
    # Find all amino groups
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    # Count the number of amino groups and return accordingly
    if len(amino_matches) >= 2:
      return True, "Contains two or more amino groups"
    else:
        return False, "Does not contain two or more amino groups"