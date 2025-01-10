"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine has a specific L-carnitine motif with a positively charged nitrogen 
    and esterified acyl group.

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
    # Looking for ester bond pattern with specific carnitine structure and L-chirality
    o_acyl_l_carnitine_pattern = Chem.MolFromSmarts("[C@@H]([C@@H](C=O)O[C@H](CN(C)(C)C)C=O)")

    if not mol.HasSubstructMatch(o_acyl_l_carnitine_pattern):
        return False, "Does not match refined O-acyl-L-carnitine pattern with L-chirality"

    return True, "Contains refined characteristic O-acyl-L-carnitine pattern with L-chirality"