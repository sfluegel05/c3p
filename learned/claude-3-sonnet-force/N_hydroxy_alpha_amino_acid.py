"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: CHEBI:26205 N-hydroxy-alpha-amino-acid
Any amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amino acid backbone pattern (C-C-N-C-C=O)
    amino_acid_pattern = Chem.MolFromSmarts("[C,$(C(-C)(-C)(-C))](-C(-N(-[#6]))(-C(-C(=O))))")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"
    
    # Look for N-hydroxy group (-N(-[OH]) or -N(-[O])(-[O]))
    n_hydroxy_pattern = Chem.MolFromSmarts("[N;X3](-[OH])[#6]")
    n_dihydroxy_pattern = Chem.MolFromSmarts("[N;X3](-[O])(-[O])[#6]")
    if not mol.HasSubstructMatch(n_hydroxy_pattern) and not mol.HasSubstructMatch(n_dihydroxy_pattern):
        return False, "No N-hydroxy group found"
    
    return True, "Contains N-hydroxy group in alpha position of amino acid backbone"