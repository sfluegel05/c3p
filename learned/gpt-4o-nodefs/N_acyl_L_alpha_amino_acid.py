"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid typically has an acyl group connected through the nitrogen
    to the backbone of an L-alpha-amino acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acyl-L-alpha-amino acid
    n_acyl_l_alpha_amino_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](C)C(=O)O")
    
    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(n_acyl_l_alpha_amino_pattern):
        return True, "Contains N-acyl group attached to L-alpha amino acid backbone"
    
    return False, "Does not match N-acyl-L-alpha-amino acid structure"

# Example usage:
test_smiles = "CC(=O)NCCCC[C@H](N)C(O)=O"  # Example SMILES for N(6)-acetyl-L-lysine
print(is_N_acyl_L_alpha_amino_acid(test_smiles))