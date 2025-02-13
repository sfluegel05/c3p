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
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a peptide bond pattern: N-C(=O)C
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)C")
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"

    # Find all peptide bond matches
    matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(matches) < 10:
        return False, f"Contains {len(matches)} amino acid residues, need at least 10 for a polypeptide"

    return True, "Contains 10 or more amino acid residues; classified as a polypeptide"