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
    
    # Improved amide linkage pattern (NC(=O)C) allowing variable attachment
    # Ensuring the carbon following the amide C=O is attached to a chain
    amide_pattern = Chem.MolFromSmarts("N[C;H2][C](=O)[C]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No acyl-amide linkage found"

    # Improved sphingosine backbone pattern
    # Including a long hydrocarbon chain, hydroxyl groups, amide linkage,
    # and often a trans double bond
    sphingosine_pattern = Chem.MolFromSmarts("[C@H]([C@@H](CO)O)C=C")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone pattern found"
        
    # Check for minimum carbon chain length for typical N-acylsphingosines
    carbon_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_atom_count < 24:
        return False, "Insufficient carbon chain length for N-acylsphingosine"

    # Ensure appropriate number of hydroxyl groups in precise locations
    correct_hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](CO)O")
    if not mol.HasSubstructMatch(correct_hydroxyl_pattern):
        return False, "Required hydroxyl groups in proper configuration not found"

    # All checks passed, classify as N-acylsphingosine
    return True, "Contains sphingosine backbone with an acyl-amide linkage"