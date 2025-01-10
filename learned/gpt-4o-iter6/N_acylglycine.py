"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is characterized by an acyl group attached to the nitrogen atom
    of glycine, which typically involves a C(=O)N-C-COOH structure where the nitrogen is
    attached to an acyl group.

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

    # Broad pattern for N-acylglycine: Look for an amide (C(=O)N) attachment to an amino acid (typically C-C(=O)O)
    n_acylglycine_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[CX4]-[CX3](=O)[OX2H1]")
    
    if mol.HasSubstructMatch(n_acylglycine_pattern):
        return True, "Molecule contains an N-acylglycine structure"

    return False, "No N-acylglycine structure found in the molecule"