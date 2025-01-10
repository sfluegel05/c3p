"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine generally has a catechol (benzene-1,2-diol) structure
    and an aminoethyl chain with possible minor variations.

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
    
    # Tighten catechol pattern: benzene ring with diol at positions 1 and 2
    catechol_pattern = Chem.MolFromSmarts("c1ccc(O)c(O)c1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No strict catechol (benzene-1,2-diol) moiety found"
    
    # Improved pattern for aminoethyl chain: flexible but specific enough for catecholamines
    aminoethyl_pattern = Chem.MolFromSmarts("NCC(O)")
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No specific aminoethyl chain found"
    
    return True, "Contains catechol structure with well-defined aminoethyl chain, classified as catecholamine"