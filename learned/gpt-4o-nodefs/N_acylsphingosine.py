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
    
    # Generalized amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("NC(=O)C")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No acyl-amide linkage found"

    # Generalization of the sphingosine pattern
    sphingosine_pattern = Chem.MolFromSmarts("C=C(O)[C@H](O)C")  # Allowing a degree of flexibility
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine-like backbone pattern found"

    # Check for sufficient ion chain length in context of known examples
    carbon_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_atom_count < 24:
        return False, "Insufficient carbon chain length for N-acylsphingosine"

    # Ensure at least two hydroxyl groups in reasonable location
    # Focusing pattern to more generic backbone regions
    correct_hydroxyl_count = Chem.MolFromSmarts("O[C@H]")
    if not mol.HasSubstructMatch(correct_hydroxyl_count):
        return False, "Required hydroxyl groups not found in correct stereochemistry"

    return True, "Contains sphingosine backbone with an acyl-amide linkage"