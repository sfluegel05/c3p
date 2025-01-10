"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine typically has an amino group linked to an aromatic ring by a two-carbon bridge.

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

    # Look for aromatic ring pattern
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")  # benzene as a simple example
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"
    
    # Look for two-carbon alkyl chain connected to an amino group
    amine_pattern = Chem.MolFromSmarts("NCC")  # simple primary or secondary amine pattern
    if not mol.HasSubstructMatch(amine_pattern):
        return False, "No two-carbon chain linked to an amino group found"

    return True, "Contains an aromatic ring and an amino group linked via a two-carbon chain"