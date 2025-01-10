"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    
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

    # Pattern for O-acyl-L-carnitine: Ester linkage with L-carnitine characteristics
    o_acyl_l_carnitine_pattern = Chem.MolFromSmarts("O[C@H](CC(=O)[O-])C[N+](C)(C)C")
    if not mol.HasSubstructMatch(o_acyl_l_carnitine_pattern):
        return False, "No O-acyl-L-carnitine pattern found"
    
    return True, "Contains O-acyl-L-carnitine pattern"