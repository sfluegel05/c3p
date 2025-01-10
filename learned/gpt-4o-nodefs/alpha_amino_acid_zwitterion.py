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
    
    # Correct pattern for detecting an alpha-amino acid zwitterion
    # An alpha carbon (C) connected to an ammonium group ([NH3+])
    # and a carboxylate group ([C(=O)[O-]])
    zwitterion_pattern = Chem.MolFromSmarts("[C@@H]([NH3+])[C](=O)[O-]")
    
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "No zwitterionic alpha carbon with NH3+ and COO- found"

    return True, "Contains alpha carbon with NH3+ and COO-, indicating alpha-amino-acid zwitterion"