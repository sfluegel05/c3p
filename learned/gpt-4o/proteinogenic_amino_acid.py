"""
Classifies: CHEBI:83813 proteinogenic amino acid
"""
from rdkit import Chem

def is_proteinogenic_amino_acid(smiles: str):
    """
    Determines if a molecule is a proteinogenic amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a proteinogenic amino acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES to create RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycine special case: it's non-chiral
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)O")
    if mol.HasSubstructMatch(glycine_pattern):
        return True, "Glycine detected, which is an achiral proteinogenic amino acid"

    # Check for standalone amino group and carboxyl group on the same alpha carbon
    alpha_amino_carboxyl_pattern = Chem.MolFromSmarts("N[C@H](C(=O)O)C")
    if not mol.HasSubstructMatch(alpha_amino_carboxyl_pattern):
        return False, "No standalone amino acid detected"

    # Check for peptide bond patterns to exclude peptides
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Peptide bond detected, molecule may be part of a peptide"

    return True, "The molecule fits the profile of a standalone proteinogenic amino acid"