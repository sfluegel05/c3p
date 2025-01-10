"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid has an acyl group attached to the nitrogen of an L-alpha-amino acid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Create a molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more robust SMARTS pattern for detecting N-acyl-donation on an amino acid structure
    # Allow for detection of acyl attachments both directly to the alpha nitrogens
    pattern_main = Chem.MolFromSmarts("C(=O)[N;X3][C@H]([C;R0])[C;R0](=O)[O;X1,O-]")
    
    # Also consider variations that appear in the side-chains of amino acids like lysine
    pattern_side_chain = Chem.MolFromSmarts("C(=O)[N;X3][a,c][C;H1,H2]")
    
    if mol.HasSubstructMatch(pattern_main):
        return True, "Matches N-acyl-L-alpha-amino acid structure (main chain)"
    
    if mol.HasSubstructMatch(pattern_side_chain):
        return True, "Matches N-acyl-L-alpha-amino acid structure (side chain)"

    return False, "Does not match N-acyl-L-alpha-amino acid structure"

# Example usage:
test_smiles = "CC(=O)NCCCC[C@H](N)C(O)=O"  # Example SMILES for N(6)-acetyl-L-lysine
print(is_N_acyl_L_alpha_amino_acid(test_smiles))