"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid has an aromatic ring and an amino acid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more comprehensive SMARTS pattern for any aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("a")
    has_aromatic_ring = mol.HasSubstructMatch(aromatic_ring_pattern)
    
    # Refine the SMARTS pattern for amino acid backbone to allow more flexibility
    # Covering [NX3] (amino group) and [CX3](=O)[OX1H] (carboxylic acid) attached to an alpha carbon
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4]([*])[CX3](=O)[OX2H1]")
    has_amino_acid_structure = mol.HasSubstructMatch(amino_acid_pattern)
    
    if has_aromatic_ring and has_amino_acid_structure:
        return True, "Contains both an aromatic ring and an amino acid backbone"
    elif not has_aromatic_ring:
        return False, "Missing aromatic ring"
    elif not has_amino_acid_structure:
        return False, "Missing amino acid structure"

    return False, "Does not match criteria for an aromatic amino acid"