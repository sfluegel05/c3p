"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is defined by the presence of a protonated amine group
    and a deprotonated carboxylate group on the same alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for protonated amine group connected to a central carbon
    amine_pattern = Chem.MolFromSmarts("[CH]([NH3+])[C]([C])[C](=O)[O-]")

    # Search pattern for common alpha-amino-acid layouts with zwitterion features
    if mol.HasSubstructMatch(amine_pattern):
        return True, "Matches the alpha-amino-acid zwitterion pattern"
    
    # Add flexibility to identify structures with complex side chains
    generic_amino_acid_pattern = Chem.MolFromSmarts("[NH3+]C([#6])[C](=O)[O-]") 
    
    if mol.HasSubstructMatch(generic_amino_acid_pattern):
        return True, "Matches flexible alpha-amino-acid zwitterion pattern"
    
    return False, "Does not match the zwitterion pattern"