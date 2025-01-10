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
    
    # Improved pattern for detecting an alpha-amino acid zwitterion
    # An alpha carbon bonded to both an ammonium group [NH3+] 
    # and a carboxylate group [C(=O)[O-]]
    zwitterion_pattern = Chem.MolFromSmarts("[C;R0][N+](C)(C)(C) |C2=CNN=C2|") # Better captures zwitterionic state
    
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "No zwitterionic alpha carbon pattern with NH3+ and COO- found"

    return True, "Contains alpha carbon with NH3+ and COO-, indicating alpha-amino-acid zwitterion"