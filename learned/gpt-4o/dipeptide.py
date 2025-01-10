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

    # Improved pattern for identifying an amide bond [N-C(=O)] with potential variations
    amide_pattern = Chem.MolFromSmarts("[$([NX3][CX3]=[O])]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)

    # Check for at least two amide bonds
    if len(amide_matches) < 2:
        return False, f"Found {len(amide_matches)} amide bonds, need at least 2 for a dipeptide"

    # Improved pattern for identifying amino acid residues
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H,$(N-C-C(=O))]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # Check for at least two amino acid residues
    if len(amino_acid_matches) < 2:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, need at least 2"

    # Check overall connectivity for typical dipeptide configuration
    # This assumes two amino acids connected through a peptide bond
    peptide_bond_count = sum(1 for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 1.33)
    if peptide_bond_count < 2:
        return False, "Insufficient connectivity pattern for dipeptide"

    return True, "Contains at least two amino-acid residues with peptide linkages (amide bonds)"