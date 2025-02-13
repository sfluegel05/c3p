"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide contains two amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Search for amide bonds (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("N-C(=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    if len(amide_matches) != 2:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 2 for a dipeptide"

    # Look for amino acid residues (N-C-C(=O))
    amino_acid_pattern = Chem.MolFromSmarts("N-C-C(=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    if len(amino_acid_matches) != 2:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, need exactly 2"

    return True, "Contains two amino-acid residues with peptide linkages"