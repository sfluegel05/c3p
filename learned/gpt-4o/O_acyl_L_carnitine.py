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
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for O-acyl-L-carnitine considering acyl group and chiral center.
    # L-Carnitine: Look for [_:1]C([_:2])([*:3])C[N+](C)(C)C with (C-C-C[N+]) configuration at the chiral center
    l_carnitine_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CC(=O)[O-])C[N+](C)(C)C")
    l_carnitine_inverse_pattern = Chem.MolFromSmarts("C(=O)O[C@@H](CC(=O)[O-])C[N+](C)(C)C")
    
    # Check for presence of the L-carnitine configuration
    if mol.HasSubstructMatch(l_carnitine_pattern):
        return True, "Correctly matches O-acyl-L-carnitine pattern with L-configuration"
    
    # Check for alternative scenario with 'C@@H' (inverse stereochemistry that might still be L)
    if mol.HasSubstructMatch(l_carnitine_inverse_pattern):
        return True, "Correctly matches O-acyl-L-carnitine pattern with inverse stereochemistry"
    
    return False, "Does not match O-acyl-L-carnitine pattern"