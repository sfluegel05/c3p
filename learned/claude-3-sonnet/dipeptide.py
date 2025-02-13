"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: CHEBI:15554 dipeptide

A dipeptide is any molecule that contains two amino-acid residues connected by peptide linkages.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dipeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bond pattern (-C(=O)-N-C-C(=O)-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)NC(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Count number of peptide bonds
    n_peptide_bonds = len(peptide_bond_matches)
    if n_peptide_bonds != 1:
        return False, f"Found {n_peptide_bonds} peptide bonds, expected 1 for a dipeptide"

    # Look for amino acid residues (N-C-C)
    amino_acid_pattern = Chem.MolFromSmarts("[N;X3;H2,H1]([C;X4])[C;X4]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # Count number of amino acid residues
    n_amino_acids = len(amino_acid_matches)
    if n_amino_acids != 2:
        return False, f"Found {n_amino_acids} amino acid residues, expected 2 for a dipeptide"

    return True, "Contains two amino acid residues connected by a peptide linkage"