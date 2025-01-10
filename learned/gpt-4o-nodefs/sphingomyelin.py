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

    # Revised phosphate-choline headgroup pattern to accommodate more variability
    phosphate_choline_pattern = Chem.MolFromSmarts("C[P;R](=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphate_choline_pattern):
        return False, "No phosphate-choline group found"

    # More general sphingoid base pattern recognizing variability
    sphingoid_base_pattern = Chem.MolFromSmarts("[C@@H](N)[C@@H](O)CC")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base found"

    # Long hydrocarbon chain pattern with allowance for unsaturation
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C(=O)CCCCCCCCCCCCCCCC")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficiently long hydrocarbon chain found"

    return True, "Structure matches sphingomyelin characteristics"