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
    
    # Zwitterion pattern for alpha-amino acids:
    # Positive nitrogen attached to one alpha carbon, the same carbon attached to a carboxylate
    # Pattern for a positively charged nitrogen attached to an alpha carbon
    pos_nitrogen_pattern = Chem.MolFromSmarts("[NH1+,NH2+,NH3+]")  # Consideration for different positive nitrogen

    # Pattern for carboxylate attached to the alpha carbon
    carboxylate_pattern = Chem.MolFromSmarts("[C](=O)[O-]")

    pos_nitrogen_matches = mol.HasSubstructMatch(pos_nitrogen_pattern)
    carboxylate_matches = mol.HasSubstructMatch(carboxylate_pattern)

    if pos_nitrogen_matches and carboxylate_matches:
        return True, "Contains alpha carbon with NH3+ and COO-, indicating alpha-amino-acid zwitterion"
    else:
        return False, "No zwitterionic alpha carbon with NH3+ and COO- found"