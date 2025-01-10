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

    # Define the phosphate-choline headgroup pattern
    phosphate_choline_pattern = Chem.MolFromSmarts("COP([O-])(=O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphate_choline_pattern):
        return False, "No phosphate-choline group found"

    # Define the sphingoid base pattern, focusing on common chiral centers and base
    sphingoid_base_pattern = Chem.MolFromSmarts("[C@@H](NC(=O)C)[C@@H](O)C[C;X4]")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base found"

    # Define a flexible long hydrocarbon chain pattern (allow for variations, including double bonds)
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("C[C@@H](O)CC([C;!H0])")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No sufficiently long hydrocarbon chain found"

    return True, "Structure matches sphingomyelin characteristics"