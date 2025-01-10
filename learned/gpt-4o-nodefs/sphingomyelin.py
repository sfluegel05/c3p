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

    # Broad phosphate-choline pattern, capturing possible variations
    phosphate_choline_pattern = Chem.MolFromSmarts("O=P(OCC[N+](C)(C)C)(O)O")
    if not mol.HasSubstructMatch(phosphate_choline_pattern):
        return False, "No phosphate-choline group found"

    # Broadened sphingoid base detection pattern
    # This pattern attempts to capture common structural variations
    sphingoid_base_pattern = Chem.MolFromSmarts("[C@@H](NC(=O)[C;X4])[C@@H](O)[C;X4]")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base found"

    # Flexible long hydrocarbon chain pattern (allow for variations)
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C(CCCCCCCCCCC)*)")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficiently long hydrocarbon chain found"

    return True, "Structure matches sphingomyelin characteristics"