"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: CHEBI:32988 oligopeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is a peptide containing a relatively small number of amino acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide bonds (-C(=O)-N-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)NC")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)

    # If no peptide bonds, it's not a peptide
    if not peptide_bonds:
        return False, "Not a peptide"

    # Count amino acid residues
    amino_acid_pattern = Chem.MolFromSmarts("[NX3H2,NX4H3+][C@H](N)C(=O)"
                                              "|"
                                              "[NX3H0,NX4H+0]([C@@H](C(=O))"
                                              "[C@@H](N)C(=O))[C@@H](N)C(=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    n_residues = len(amino_acid_matches)

    # Oligopeptides typically have < 10 amino acids
    if n_residues < 10:
        return True, f"Contains {n_residues} amino acid residues"
    else:
        return False, f"Contains {n_residues} amino acid residues, likely not an oligopeptide"