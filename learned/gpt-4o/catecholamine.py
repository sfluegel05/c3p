"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine generally has a catechol (benzene-1,2-diol) structure
    and a 2-aminoethyl chain with possible variations.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a flexible catechol pattern (benzene-1,2-diol) allowing for substituted forms
    flexible_catechol_pattern = Chem.MolFromSmarts("c1c(O)c(O)c(*)c(*)c1")
    if not mol.HasSubstructMatch(flexible_catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) moiety found"
    
    # Look for an aminoethyl-like pattern (-[CH2]-[CH2]-[N]-)
    # Flexible enough to include variations such as N-methyl or additional alkyl groups
    aminoethyl_variation_pattern = Chem.MolFromSmarts("[NX3][CH2][CH2]")
    if not mol.HasSubstructMatch(aminoethyl_variation_pattern):
        return False, "No flexible 2-aminoethyl-like chain found"
    
    return True, "Contains catechol structure with 2-aminoethyl-like chain, classified as catecholamine"