"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: CHEBI:63850 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is defined as an aldehyde with a long carbon chain, typically
    derived from the reduction of a fatty acid.

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

    # Look for aldehyde group pattern (C=O at terminal position)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No terminal aldehyde group found"
    
    # Look for a sufficiently long carbon chain with flexibility for double/triple bonds and branching
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        # If strict chain pattern is not found, check for unsaturated or branched alternatives
        unsaturated_chain_pattern = Chem.MolFromSmarts("[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]~[C,c]")
        if not mol.HasSubstructMatch(unsaturated_chain_pattern):
            return False, "Does not have a sufficient long carbon chain"

    # Count total number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Adjust carbon count threshold considering smaller examples and flexibility
    if c_count < 6:
        return False, f"Too few carbons for a typical fatty aldehyde (found {c_count})"

    return True, "Contains terminal aldehyde group with a suitable long carbon chain"