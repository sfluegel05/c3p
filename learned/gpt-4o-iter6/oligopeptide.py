"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is defined as a peptide containing a relatively small number of amino acids (typically 2-20).

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
    
    # Check for presence of peptide bonds (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    num_peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))
    if num_peptide_bonds < 1:
        return False, "No peptide bonds identified"
    
    # Determine if the number of peptide bonds falls within oligopeptide range
    if num_peptide_bonds > 19:
        return False, f"Too many peptide bonds for an oligopeptide: {num_peptide_bonds}"
    
    # Assess for alpha amino acid patterns (common in peptides)
    alpha_amino_pattern = Chem.MolFromSmarts("[C@@H](N)C(=O)O")
    num_alpha_amino_residues = len(mol.GetSubstructMatches(alpha_amino_pattern))
    if num_alpha_amino_residues < 2:
        return False, "Too few alpha amino acid residues for an oligopeptide"

    return True, "Oligopeptide structure identified with appropriate number of peptide bonds"