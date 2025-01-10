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

    # Look for the glycine moiety, accounting for stereochemistry.
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine moiety found"

    # Look for the acyl group directly attached to the nitrogen of the glycine moiety.
    # Allows variable carbon attachment to cover more diverse N-acyl variations.
    n_acyl_pattern = Chem.MolFromSmarts("N[C;R0][CX4,CX3](=O)")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No correctly attached acyl group found"

    return True, "Contains glycine moiety linked to an acyl group with proper attachment"

# Example usage:
# test_smiles = "CCCCCCCCCCCC(=O)NCC(O)=O"  # SMILES for N-dodecanoylglycine
# is_N_acylglycine(test_smiles)