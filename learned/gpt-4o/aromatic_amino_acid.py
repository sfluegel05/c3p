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
    
    # SMARTS pattern for aromatic ring (e.g., benzene)
    aromatic_ring_pattern = Chem.MolFromSmarts("a1aaaaa1")
    has_aromatic_ring = mol.HasSubstructMatch(aromatic_ring_pattern)
    
    # SMARTS pattern for amino acid backbone (NH2-CH(R)-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H][CX3](=O)[OX2H1]")
    has_amino_acid_structure = mol.HasSubstructMatch(amino_acid_pattern)
    
    if has_aromatic_ring and has_amino_acid_structure:
        return True, "Contains both an aromatic ring and an amino acid backbone"
    elif not has_aromatic_ring:
        return False, "Missing aromatic ring"
    elif not has_amino_acid_structure:
        return False, "Missing amino acid structure"

    return False, "Does not match criteria for an aromatic amino acid"