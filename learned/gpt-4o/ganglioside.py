"""
Classifies: CHEBI:28892 ganglioside
"""
from rdkit import Chem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a molecule composed of a glycosphingolipid (ceramide and oligosaccharide) 
    with one or more sialic acids linked on the sugar chain.

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

    # Define updated patterns
    ceramide_pattern = Chem.MolFromSmarts("NC([C@@H](CO)O)C(=O)[C@H](CCCCCCCCCCCCCCCC)O")  # Sphingosine backbone and fatty acid
    # Simplified oligosaccharide chain with common sugars (could be adjusted for specific sugars)
    oligosaccharide_pattern = Chem.MolFromSmarts("O[C@H]([C@@H](O)CO)O")  # Glc/Gal building block
    # Sialic acid pattern to detect presence - allowing for potential N-acetyl group
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@@H](O)[C@H](C(=O)O)[C@@H](O)C")  # Expanded sialic acid matcher

    # Check for ceramide presence
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Check for an oligosaccharide chain
    sugar_matches = mol.GetSubstructMatches(oligosaccharide_pattern)
    if len(sugar_matches) < 2:  # Assume minimum 2 sugar units
        return False, f"Insufficient sugar moieties, found {len(sugar_matches)}"

    # Check for at least one sialic acid unit
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_matches) == 0:
        return False, "No sialic acid found"

    return True, "Contains ceramide backbone with glycosphingolipid and sialic acid(s)"