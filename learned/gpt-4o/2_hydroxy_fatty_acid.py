"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy functional group in the alpha- or 2-position
    relative to the carboxylic acid group, with a sufficiently long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns
    # Main components of 2-hydroxy fatty acid
    alpha_hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("O[C@@H](C)C(=O)O")
    # A simpler non-stereospecific alpha hydroxy group check
    non_stereo_alpha_hydroxy_pattern = Chem.MolFromSmarts("OCC(=O)O")
    
    # Minimum number of carbons in the chain (e.g., 6 for hexanoic acids)
    # [CH2]6 for a chain of at least 6 methylene (-CH2-) groups
    min_chain_length_pattern = Chem.MolFromSmarts("[C]([C])([C])[C]")

    # Check if it has alpha hydroxy fatty acid pattern
    if not mol.HasSubstructMatch(alpha_hydroxy_fatty_acid_pattern) and not mol.HasSubstructMatch(non_stereo_alpha_hydroxy_pattern):
        return False, "No hydroxy group adjacent to the carboxylic acid group"

    # Check if it has a minimal chain length
    if not mol.HasSubstructMatch(min_chain_length_pattern):
        return False, "Insufficient carbon chain length typical for a fatty acid"

    return True, "Contains hydroxy group at 2-position with appropriate aliphatic chain and carboxylic acid group"

# Examples of correct SMILES for debugging
# is_2_hydroxy_fatty_acid("CCCCCCCC(O)C(O)=O")  # Should return True
# is_2_hydroxy_fatty_acid("OC(C(O)=O)C")  # Should return False because it is too short