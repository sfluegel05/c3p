"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid contains a 3beta-hydroxy group and a double bond between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns
    ## More precise pattern for 3beta-hydroxy
    hydroxy_beta_pattern = Chem.MolFromSmarts("[C@@H]([O])[C;R1]")
    
    ## Properly specify Delta(5) double bond in steroid context
    # it usually involves specific connectivity
    delta5_specific_pattern = Chem.MolFromSmarts("C1[C@@]2(C=C[C@@H]3)C[C@@H](C)(CCC3)CC2C(C1)O")  # generalization
    
    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_beta_pattern):
        return False, "No 3beta-hydroxy group found"
    
    # Check for Delta(5) double bond match
    if not mol.HasSubstructMatch(delta5_specific_pattern):
        return False, "No Delta(5) double bond found (C5-C6 double bond)"

    return True, "Contains 3beta-hydroxy group with Delta(5) double bond (C5-C6 double bond)"