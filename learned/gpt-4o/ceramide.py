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

    # Sphingoid base: long-chain amino alcohol with a generalized pattern
    sphingoid_pattern = Chem.MolFromSmarts("N[C@H](CO)[C@@H](O)CC")
    
    # Amide linkage pattern - simplified
    amide_pattern = Chem.MolFromSmarts("N[C@H](O)C(=O)C")  # More general, detect amide

    # Check for sphingoid base
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "Sphingoid base pattern not found"

    # Check for amide linkage
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Amide linkage pattern not found"
    
    # Ensure carbon chain within reasonable range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14 or c_count > 26:
        return False, f"Number of carbon atoms ({c_count}) is outside typical ceramide range (14-26)"
    
    # Check for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")  # Recognize -OH
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No common hydroxyl groups found"

    return True, "Contains characteristic sphingoid base with an amide-linked fatty acid chain"