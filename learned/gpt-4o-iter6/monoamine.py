"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is an arylamino compound containing one amino group connected
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
        return False, "Invalid SMILES string"

    # SMARTS pattern for amino group connected by a two-carbon chain to an aromatic carbon
    monoamine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;R0]-C-C-[Ar]")  # N represented as NH2 or NH

    # Check for a two-carbon chain with amino group connected to an aromatic system
    if mol.HasSubstructMatch(monoamine_pattern):
        return True, "Contains an amino group connected to aromatic ring by a two-carbon chain"

    return False, "Missing correct structure configuration for monoamine classification"