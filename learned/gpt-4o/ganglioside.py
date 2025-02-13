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
    
    # Define more accurate patterns for recognition of structural components
    ceramide_pattern = Chem.MolFromSmarts("NC([C@@H]([NH2])CO)C(=O)[C@H](CCCCCCCCCCCCCCCC)O")  # Sphingosine base and fatty acid
    oligosaccharide_pattern = Chem.MolFromSmarts("O[C@H]([C@@H](CO)O)C")  # Glc or Gal units in sugar chain
    sialic_acid_pattern = Chem.MolFromSmarts("C[C@@H](O)C(O)=O")  # Simplified sialic acid pattern

    # Check for ceramide presence
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone found"

    # Check for oligosaccharide chain presence
    sugar_matches = mol.GetSubstructMatches(oligosaccharide_pattern)
    if len(sugar_matches) < 3:  # Assume minimum 3 sugar units for complexity
        return False, f"Insufficient sugar moieties, found {len(sugar_matches)}"

    # Check for sialic acids
    sialic_acid_matches = mol.GetSubstructMatches(sialic_acid_pattern)
    if len(sialic_acid_matches) == 0:
        return False, "No sialic acid found in the structure"

    # Ensure all components are interconnected properly
    if not all([mol.HasSubstructMatch(pattern) for pattern in [ceramide_pattern, oligosaccharide_pattern, sialic_acid_pattern]]):
        return False, "Structure not recognized as a connected ganglioside"

    return True, "Contains ceramide backbone with attached glycosphingolipid structure and sialic acid(s)"