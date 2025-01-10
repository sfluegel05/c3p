"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible peptide bond pattern considering standard amide linkage:
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N[C@@H,C@H]")  # Broadens the stereochemistry pattern
    if peptide_bond_pattern is None:
        return None, "Failed to parse peptide bond SMARTS pattern"

    # Check for recognizable peptide bonds in the molecule
    peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    amino_acid_count = len(peptide_matches)

    # Assume bonds count relates to amino acids, since patterns can repeat
    if amino_acid_count < 10:
        return False, f"Contains {amino_acid_count} peptide bonds, need at least 10 for a polypeptide"

    # Check for terminal groups (N and C-terminals) to ensure complete polypeptide structure
    # Assume that we need NH2 (or similar) and COOH (or similar), encapsulated in a pattern search:
    n_terminal_pattern = Chem.MolFromSmarts("[NH2,NH][CH,C][=O]")
    c_terminal_pattern = Chem.MolFromSmarts("C(=O)[OH]")

    if not (mol.HasSubstructMatch(n_terminal_pattern) and mol.HasSubstructMatch(c_terminal_pattern)):
        return False, "Does not contain typical terminal groups expected in polypeptides"

    return True, "Contains 10 or more peptide bonds; classified as a polypeptide with proper terminals"