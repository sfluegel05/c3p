"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: CHEBI:17478 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde with a long carbon chain, typically arising from the reduction of a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an aldehyde group (C=O with at least one hydrogen)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#1]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        # Try a more general pattern for aldehydes
        aldehyde_pattern = Chem.MolFromSmarts("[CX3](=O)[#1]")
        if not mol.HasSubstructMatch(aldehyde_pattern):
            return False, "No aldehyde group found"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Carbon chain too short to be a fatty aldehyde"

    # Check for long carbon chain (fatty aldehyde characteristic)
    # A fatty aldehyde typically has a carbon chain length of 6 or more
    if c_count >= 6:
        return True, "Contains a long carbon chain with an aldehyde group"
    else:
        return False, "Carbon chain too short to be a fatty aldehyde"