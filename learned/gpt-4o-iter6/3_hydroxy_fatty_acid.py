"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a hydroxy functional group at the beta- or 3-position
    from the carboxylic acid and is characterized by a long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Locate the carboxylic acid group using a SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for a hydroxy group at the 3-position (beta-position) from the carboxylic acid
    hydroxy_at_3_pattern = Chem.MolFromSmarts("C(=O)[OH][CH2][CH2][O]")
    if not mol.HasSubstructMatch(hydroxy_at_3_pattern):
        return False, "No hydroxy group at the 3-position"

    # Determine length of carbon chain - Fatty acids typically have chains with 8 or more carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 8:
        return False, "Carbon chain too short for a fatty acid"

    return True, "Contains a hydroxy group at the 3-position and a sufficient carbon chain"