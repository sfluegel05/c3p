"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid should have an aromatic ring and a simple amino acid backbone.

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
    
    # Detect any aromatic ring within the molecule
    aromatic_ring_pattern = Chem.MolFromSmarts("a")
    has_aromatic_ring = mol.HasSubstructMatch(aromatic_ring_pattern)
    
    # Use a SMARTS pattern for a simple amino acid backbone:
    # - [NX3;H2] for a primary amine group
    # - [CX4;H] for the alpha carbon with a side chain
    # - [CX3](=O)[OX1H] representing the carboxylic acid group
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2][CX4H1]([*])[CX3](=O)[OX1H]")
    has_amino_acid_backbone = mol.HasSubstructMatch(amino_acid_pattern)
    
    # Ensure there's only one primary amine and one carboxylic acid without intervening peptide bonds
    is_single_amino_acid = len(mol.GetSubstructMatches(amino_acid_pattern)) == 1

    if has_aromatic_ring and has_amino_acid_backbone and is_single_amino_acid:
        return True, "Molecule is an aromatic amino acid with one aromatic ring and one simple amino acid backbone"
    elif not has_aromatic_ring:
        return False, "Missing aromatic ring"
    elif not has_amino_acid_backbone:
        return False, "Missing amino acid structure or contains peptide bonds"

    return False, "Does not match criteria for an aromatic amino acid"