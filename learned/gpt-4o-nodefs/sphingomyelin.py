"""
Classifies: CHEBI:64583 sphingomyelin
"""
from rdkit import Chem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate-choline pattern
    phosphate_choline_pattern = Chem.MolFromSmarts("O=P(OCC[N+](C)(C)C)(O)O")
    if not mol.HasSubstructMatch(phosphate_choline_pattern):
        return False, "No phosphate-choline group found"

    # Look for sphingoid base (long chain amino alcohol)
    sphingoid_base_pattern = Chem.MolFromSmarts("[C@@H](NC(=O)C)[C@@H](O)C(CCCCCCCCCCCCC)")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base found"

    # Verify the presence of a long hydrocarbon chain
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCCC)")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No long hydrocarbon chain found"

    return True, "Structure matches sphingomyelin characteristics"