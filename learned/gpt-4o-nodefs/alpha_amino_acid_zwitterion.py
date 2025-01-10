"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
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

    # Check for presence of an alpha carbon connected to NH3+ and COO-
    alpha_carbon_pattern = Chem.MolFromSmarts("[C@H]([NH3+])[C](=[O])[O-]")
    alpha_carbon_alt_pattern = Chem.MolFromSmarts("[C@]([NH3+])[C](=[O])[O-]") # To handle relative chirality

    if not (mol.HasSubstructMatch(alpha_carbon_pattern) or mol.HasSubstructMatch(alpha_carbon_alt_pattern)):
        return False, "No alpha carbon with NH3+ and COO- found"

    # If reached here, the molecule contains the necessary functional groups
    return True, "Contains alpha carbon with NH3+ and COO-, indicating alpha-amino-acid zwitterion"