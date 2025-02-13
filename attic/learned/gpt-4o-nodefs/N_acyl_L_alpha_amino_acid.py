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
    
    # Define the SMARTS pattern for an N-acyl-L-alpha-amino acid
    # Pattern: Acyl group (O=C-) bonded to Nitrogen (N) bonded to chiral alpha carbon [C@H]
    pattern = Chem.MolFromSmarts("N([C@H])C(=O)")
    
    # Check if the molecule matches the N-acyl-L-alpha-amino acid pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Matches N-acyl-L-alpha-amino acid structure"
    else:
        return False, "Does not match N-acyl-L-alpha-amino acid structure"

# Example usage:
smiles_example = "C[C@H](NC(=O)Cc1c[nH]c2ccccc12)C(O)=O"  # Example SMILES
print(is_N_acyl_L_alpha_amino_acid(smiles_example))