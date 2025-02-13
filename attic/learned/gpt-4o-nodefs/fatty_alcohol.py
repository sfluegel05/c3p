"""
Classifies: CHEBI:24026 fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol typically has a long carbon chain with one or more hydroxyl (-OH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one hydroxyl group (-OH)
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No hydroxyl group found"

    # Check for long carbon chain - at least 6 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6:
        return False, "Carbon chain too short"

    return True, "Contains long carbon chain with hydroxyl group(s)"