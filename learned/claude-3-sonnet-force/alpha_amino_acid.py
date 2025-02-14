"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:33709 alpha-amino acid

An alpha-amino acid is defined as an amino acid in which the amino group is 
located on the carbon atom at the position alpha to the carboxy group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define alpha-amino acid pattern
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N;H2]([C;H2,H1])([C;H2,H1])([C;H2,H1])[C;H1]([N;H2,H1])([C;H1])(=O)[O;H1]")
    
    # Check if molecule matches the pattern
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    
    # Check for zwitterionic forms
    zwitterion_pattern = Chem.MolFromSmarts("[N+;H3]([C;H2,H1])([C;H2,H1])([C;H2,H1])[C;H1]([N+;H2])([C;H1])(=[O-1])")
    zwitterion_matches = mol.GetSubstructMatches(zwitterion_pattern)
    
    if matches or zwitterion_matches:
        return True, "Contains an alpha-amino acid substructure"
    
    return False, "No alpha-amino acid substructure found"