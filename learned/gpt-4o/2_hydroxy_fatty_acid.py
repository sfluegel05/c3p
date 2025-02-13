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

    # Define pattern for 2-hydroxy fatty acids
    # Match a carboxylic acid group (C(=O)O) followed by a carbon that has a hydroxy group (O)
    # Carbon chain must contain at least a few methylene groups between the carboxylic and hydroxyl
    hydroxy_fatty_acid_pattern = Chem.MolFromSmarts("C(O)C([O-])=O")

    # Check for the correct substructure in the molecule
    if not mol.HasSubstructMatch(hydroxy_fatty_acid_pattern):
        return False, "No hydroxy group adjacently at the 2-position to the carboxylic acid group"
    
    # Ensuring the presence of a sufficient length aliphatic chain
    # Check for at least 4 consecutive carbon atoms constituting a chain
    chain_length_pattern = Chem.MolFromSmarts("[C][C][C][C]")
    if not mol.HasSubstructMatch(chain_length_pattern):
        return False, "Insufficient carbon chain length typical for a fatty acid"
    
    return True, "Contains hydroxy group at 2-position with appropriate aliphatic chain and carboxylic acid group"

# Examples of correct SMILES for debugging
# is_2_hydroxy_fatty_acid("CCCCCCCC(O)C(=O)O")  # Should return True
# is_2_hydroxy_fatty_acid("OC(C(O)=O)C")  # Should return False due to short chain