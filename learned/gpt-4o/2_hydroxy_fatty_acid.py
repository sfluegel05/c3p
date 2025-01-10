"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: CHEBI:17855 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy functional group in the alpha- or 2-position relative to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    
    # Check for carboxylic acid group
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Define a pattern for 2-hydroxy group (alpha position to carboxylic acid)
    # [O;H1] is a single hydroxy group, attach to carbon that is alpha to carboxylic acid group
    hydroxy_alpha_pattern = Chem.MolFromSmarts("C([O;H1])C(=O)O")

    # Check for 2-hydroxy group in alpha position
    if not mol.HasSubstructMatch(hydroxy_alpha_pattern):
        return False, "No 2-hydroxy group found in the alpha position"
    
    # Verification of fatty acid-like long chain is more complex
    # We will simply check if the chain contains at least 6 carbons
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    
    if carbon_count < 6:
        return False, "Insufficient carbon chain length for fatty acid"

    return True, "Contains 2-hydroxy group in the alpha position of a fatty acid chain"