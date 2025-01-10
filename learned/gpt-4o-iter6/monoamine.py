"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an arylamino compound which contains one amino group connected
    to an aromatic ring by a two-carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Update SMARTS pattern: Detect primary or secondary amine connected to a two-carbon chain leading to an aromatic ring
    two_carbon_amino_pattern = Chem.MolFromSmarts("[NX3;R0]-[CH2]-[CH2]-[c]")  # 'c' is lowercase to ensure it's aromatic
    
    # Check for a two-carbon chain with amino group connected to an aromatic system
    if mol.HasSubstructMatch(two_carbon_amino_pattern):
        return True, "Contains an amino group connected to aromatic ring by a two-carbon chain"

    return False, "Missing correct structure configuration for monoamine classification"