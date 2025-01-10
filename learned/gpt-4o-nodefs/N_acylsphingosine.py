"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for acyl-amide linkage pattern (NC(=O)C)
    amide_pattern = Chem.MolFromSmarts("NC(=O)[C]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No acyl-amide linkage found"
    
    # Look for improved sphingosine backbone pattern
    # Pattern should include long hydrocarbon with trans double bond, hydroxyl groups and an amide linkage
    # Consider a general sphingosine pattern
    sphingosine_pattern = Chem.MolFromSmarts("[C@H]([C@H](CO)O)C=C")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone pattern found"
        
    # Ensure there are two correct hydroxyl groups checking the relative position
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](CO)O")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Required two hydroxyl groups pattern not found"
    
    # Check for sufficient carbon chain length typical of sphingosine and acyl chain
    carbon_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_atom_count < 24:
        return False, "Insufficient carbon chain length for N-acylsphingosine"

    # If all criteria meet, return True
    return True, "Contains sphingosine backbone with an acyl-amide linkage"