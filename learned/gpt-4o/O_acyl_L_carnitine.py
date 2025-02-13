"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine is an O-acylcarnitine in which the carnitine component has L-configuration.

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
    
    # Look for L-carnitine structure with ester linkage
    # The ester pattern and L-carnitine configuration pattern ([C@H] or [C@@H]) are demonstrated
    ester_pattern = Chem.MolFromSmarts("O[C@H](CC([O-])=O)C[N+](C)(C)C")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Matches O-acyl-L-carnitine pattern"
    
    # Invert configuration if '[C@@H]' pattern is used
    ester_pattern_inverse = Chem.MolFromSmarts("O[C@@H](CC([O-])=O)C[N+](C)(C)C")
    if mol.HasSubstructMatch(ester_pattern_inverse):
        return True, "Matches O-acyl-L-carnitine pattern with inverse stereochemistry"

    return False, "Does not match O-acyl-L-carnitine pattern"