"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine typically includes a catechol (benzene-1,2-diol) structure
    and a 2-aminoethyl chain.

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
    
    # Look for catechol pattern (benzene-1,2-diol)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) moiety found"
    
    # Look for 2-aminoethyl chain pattern (-CCN)
    aminoethyl_pattern = Chem.MolFromSmarts("CCN")
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No 2-aminoethyl chain found"
    
    return True, "Contains catechol structure with 2-aminoethyl chain, classified as catecholamine"