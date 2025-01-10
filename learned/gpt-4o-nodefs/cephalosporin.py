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

    # Define an improved SMARTS pattern for the cephalosporin core
    # Recognition of β-lactam ring (four-membered) fused with a six-membered dihydrothiazine ring
    cephalosporin_pattern = Chem.MolFromSmarts("C1C2C(N1C(=O)O)SCCC3=CSC=C23")

    # Check for cephalosporin core structure
    if not mol.HasSubstructMatch(cephalosporin_pattern):
        return False, "No cephalosporin-like core structure found"

    return True, "Contains cephalosporin-like core structure"