"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sphingoid base: long-chain amino alcohol with a general flexible pattern
    sphingoid_pattern = Chem.MolFromSmarts("N[C@@H](CO)C(O)")  # Include at least basic components
    
    # Amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("N[C@@H]C(=O)C")  # Capture amide linkage specifically related to sphingoid
    
    # Check for amide linkage & sphingoid base
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "Sphingoid base pattern not found"
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Amide linkage pattern not found"
    
    # Ensure sufficient long chain - count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14 or c_count > 26:
        return False, f"Number of carbon atoms ({c_count}) is outside typical ceramide range (14-26)"
    
    # Check for at least one hydroxyl group:
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Recognize -OH
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl groups found, common in ceramides"

    return True, "Contains characteristic sphingoid base with an amide-linked fatty acid chain"