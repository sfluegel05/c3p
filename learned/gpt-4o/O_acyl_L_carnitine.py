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
    
    # Define a SMARTS pattern for O-acyl-L-carnitine considering acyl group 
    # and specific L-carnitine chiral configuration with possible isotopes.
    l_carnitine_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CC(=O)[O-])C[N+](C)(C)C")
    l_carnitine_deuterated_pattern = Chem.MolFromSmarts("C(=O)O[C@H](CC(=O)[O-])C[N+](C)(C)C([2H])")
    
    # Check for the key feature of L-carnitine involving the correct stereochemistry
    if mol.HasSubstructMatch(l_carnitine_pattern):
        return True, "Correctly matches O-acyl-L-carnitine pattern with L-configuration"
    
    # Additional check for deuterated isotopes or modified configurations
    if mol.HasSubstructMatch(l_carnitine_deuterated_pattern):
        return True, "Correctly matches O-acyl-L-carnitine pattern including isotope features"

    return False, "Does not match O-acyl-L-carnitine pattern"