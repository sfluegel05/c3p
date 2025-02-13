"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Look for acyl-amide linkage pattern (NC(=O))
    amide_pattern = Chem.MolFromSmarts("NC(=O)")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No acyl-amide linkage found"
    
    # Look for sphingosine backbone which includes a long hydrocarbon chain with double bonds
    # Pattern: chain with \C=C\ and adjacent OH groups
    sphingosine_pattern = Chem.MolFromSmarts("[C@@H](O)[C@H](CO)\C=C\[*]")
    if not mol.HasSubstructMatch(sphingosine_pattern):
        return False, "No sphingosine backbone pattern found"
    
    # Validate if there are two hydroxyl groups (attached to consecutive carbons)
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H](O)")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Required two hydroxyl groups pattern not found"
    
    # Check for sufficient carbon chain length typical of sphingosine and acyl chain
    carbon_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_atom_count < 18:
        return False, "Insufficient carbon chain length for N-acylsphingosine"

    return True, "Contains sphingosine backbone with an acyl-amide linkage"