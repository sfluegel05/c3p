"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is a compound where a glycine moiety is linked to an acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the glycine pattern (NCC(=O)O)
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine moiety found"

    # Determine if the nitrogen has a carbonyl group (acyl group) attached
    n_acyl_pattern = Chem.MolFromSmarts("N-C(=O)")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No acyl group attached to the glycine nitrogen"

    return True, "Contains glycine moiety linked to an acyl group"

# Example usage of the function
# test_smiles = "CCCCCCCCCCCC(=O)NCC(O)=O"  # SMILES for N-dodecanoylglycine
# is_N_acylglycine(test_smiles)