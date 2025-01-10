"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide consists of exactly four amino acids linked by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a peptide bond SMARTS pattern; focusing on N-C(=O)-C linkage
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Counts the peptide bonds and ensures there are exactly 3
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3"

    # Check for 4 central or alpha carbons linked to amino groups
    alpha_pattern = Chem.MolFromSmarts("[NX3][CX4]")
    alpha_matches = mol.GetSubstructMatches(alpha_pattern)

    # Each alpha amino acid should contain these patterns for N-C linkage, count to see proper distribution
    if len(alpha_matches) != 4:
        return False, f"Found {len(alpha_matches)} N-C(alpha) patterns, need exactly 4"

    return True, "Contains exactly four amino acids linked by three peptide bonds"