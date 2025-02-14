"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a molecule composed of a glycosphingolipid (ceramide and oligosaccharide) with one or more sialic acids linked.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for recognition of structural components
    ceramide_pattern = Chem.MolFromSmarts("[NH]C([#6])C(=O)O[#6]")  # Generalized ceramide pattern
    sugar_pattern = Chem.MolFromSmarts("[C@H]1([O][C@@H]([C@H](O)[C@H](O)[C@@H]1O)O)O")  # Simplified sugar pattern
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@@H]1C[C@H](O)[C@@H]([NH])C(=O)O[C@@H]1O")  # General sialic acid pattern

    # Check for ceramide presence
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Check for adequate sugar chain
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 2:
        return False, f"Insufficient sugar moieties, found {len(sugar_matches)}"

    # Check for sialic acids
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_matches) == 0:
        return False, "No sialic acid found in the structure"

    # Ensure the structure is a connected ganglioside
    if not all([mol.HasSubstructMatch(pattern) for pattern in [ceramide_pattern, sugar_pattern, sialic_acid_pattern]]):
        return False, "Structure not recognized as a connected ganglioside"

    return True, "Contains ceramide backbone with attached glycosphingolipid structure and sialic acid(s)"