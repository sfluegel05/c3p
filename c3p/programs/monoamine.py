"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine contains one amino group connected to an aromatic ring by a two-carbon chain.

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

    # Pattern for amino group connected by a two-carbon chain to an aromatic ring
    monoamine_pattern = Chem.MolFromSmarts("[$(cc)]CCN")  # Extended SMARTS might be needed for enriched specificity
    
    # Find substructure matches
    if mol.HasSubstructMatch(monoamine_pattern):
        return True, "Contains amino group connected to an aromatic ring by a two-carbon chain"
    else:
        return False, "Does not contain the required monoamine structure"