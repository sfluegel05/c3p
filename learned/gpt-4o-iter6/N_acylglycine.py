"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is characterized by an acyl group attached to the nitrogen atom
    of glycine, typically involving a C(=O)N-C-C(=O)O structure where the nitrogen is
    attached directly to an acyl group.

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

    # Improved pattern for N-acylglycine: Specifically look for the N-acyl group attached to glycine
    # The pattern represents: Acyl group (R-C(=O)-) attached to nitrogen (N-C-C(=O)O)
    n_acylglycine_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4H2][CX3](=O)[O]")
    
    if mol.HasSubstructMatch(n_acylglycine_pattern):
        return True, "Molecule contains an N-acylglycine structure"

    return False, "No N-acylglycine structure found in the molecule"