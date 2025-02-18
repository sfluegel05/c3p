"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is characterized by a carboxylic acid group and one or more amino groups 
    in proximity, forming a recognizable amino acid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classifiable as an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refine carboxylic acid search with some biodegradability
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Refine amino group pattern to match typical amino acid nitrogens
    amino_group_pattern = Chem.MolFromSmarts("[NX3,NX4+][H,D1,D2]")
    amino_group_matches = mol.GetSubstructMatches(amino_group_pattern)
    if len(amino_group_matches) < 1:
        return False, "No amino groups found"

    # Validate the connection between amino and carboxylic groups to mimic a standard amino acid chain
    for match in amino_group_matches:
        atom = mol.GetAtomWithIdx(match[0])
        for neighbor in atom.GetNeighbors():
            if neighbor.HasSubstructMatch(carboxylic_acid_pattern):
                return True, "Contains carboxylic acid and amino group(s) in correct configuration"

    return False, "No amino acid backbone structure found"