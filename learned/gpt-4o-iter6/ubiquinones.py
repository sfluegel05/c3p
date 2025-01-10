"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones are benzoquinones with methoxy groups and a polyisoprenoid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the benzoquinone pattern
    benzoquinone_pattern = Chem.MolFromSmarts("O=C1C=C(C(=O)C=C1)C")
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No benzoquinone moiety found"

    # Check for methoxy groups (at least 2 for ubiquinones)
    methoxy_pattern = Chem.MolFromSmarts("COC")
    if len(mol.GetSubstructMatches(methoxy_pattern)) < 2:
        return False, "Insufficient methoxy groups found"
    
    # Check for long polyisoprenoid chain
    isoprenoid_chain_pattern = Chem.MolFromSmarts("[C](=C)[CH1]CCCC=C")  # Simplified pattern for long chain
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_chain_pattern)
    if len(isoprenoid_matches) < 1:
        return False, "No long polyisoprenoid chain at position 6"

    return True, "Contains a benzoquinone moiety with methoxy groups and a polyisoprenoid side chain"