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

    # Look for the glycine moiety (NCC(=O)O)
    glycine_pattern = Chem.MolFromSmiles("NCC(=O)O")
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine moiety found"

    # Look for the acyl group directly attached to glycine's nitrogen
    # Modify pattern to detect acyl variant: N-C-C(=O)
    n_acyl_pattern = Chem.MolFromSmarts("N-[#6X3](=O)-[#6]")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No correctly attached acyl group found"

    return True, "Contains glycine moiety linked to an acyl group with proper attachment"

# Example usage of the function
# test_smiles = "CCCCCCCCCCCC(=O)NCC(O)=O"  # SMILES for N-dodecanoylglycine
# is_N_acylglycine(test_smiles)