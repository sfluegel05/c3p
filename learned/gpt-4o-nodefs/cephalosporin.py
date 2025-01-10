"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are a class of β-lactam antibiotics featuring a β-lactam
    ring fused with a dihydrothiazine ring, often with additional functional
    groups at specified positions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cephalosporin, False otherwise
        str: Reason for classification or failure
    """
    
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern including general cephalosporin core features
    # Pattern reflects the 7-ACA (7-aminocephalosporanic acid) structure
    cephalosporin_pattern = Chem.MolFromSmarts("C1C(=C(N2C(C1S[C@H]2[C@@H](C(=O)O)N)=O)C(=O)O)")

    # Check for cephalosporin core structure
    if not mol.HasSubstructMatch(cephalosporin_pattern):
        return False, "No cephalosporin-like core structure found"

    return True, "Contains cephalosporin-like core structure"