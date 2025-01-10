"""
Classifies: CHEBI:23066 cephalosporin
"""
from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    Cephalosporins are a class of β-lactam antibiotics featuring a β-lactam
    ring fused with a dihydrothiazine ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cephalosporin, False otherwise
        str: Reason for classification or failure
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an improved SMARTS pattern for the cephalosporin core
    # This pattern includes beta-lactam fused to dihydrothiazine with stereochemistry acknowledgment
    cephalosporin_pattern = Chem.MolFromSmarts("C1C([C@@H]2SCC1)=C(N2)C(=O)O")

    # Check for cephalosporin ring structure
    if not mol.HasSubstructMatch(cephalosporin_pattern):
        return False, "No cephalosporin-like core structure found"

    return True, "Contains cephalosporin-like core structure"