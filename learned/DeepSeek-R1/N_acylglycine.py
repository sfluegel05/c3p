"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: CHEBI: N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is an N-acyl-amino acid where the amino acid is glycine.

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
    
    # Define SMARTS pattern for N-acylglycine core structure:
    # Acyl group (C=O) attached to nitrogen, which is connected to CH2-C(=O)OH (glycine part)
    pattern = Chem.MolFromSmarts('[NX3]([CX3](=O))[CH2][CX3](=O)[OX2H1]')
    
    if mol.HasSubstructMatch(pattern):
        return True, "Contains N-acylglycine core structure"
    
    return False, "Does not match N-acylglycine pattern"