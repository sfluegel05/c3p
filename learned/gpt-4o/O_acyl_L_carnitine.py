"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine has an L-carnitine core esterified with carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to detect O-acyl-L-carnitine
    pattern = Chem.MolFromSmarts("O=C([O;D2][C@H](C[N+](C)(C)C)CC(=O)[O-])")
    
    if not mol.HasSubstructMatch(pattern):
        return False, "Does not match O-acyl-L-carnitine pattern with L-chirality"

    return True, "Contains characteristic O-acyl-L-carnitine pattern with L-chirality"