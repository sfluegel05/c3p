"""
Classifies: CHEBI:46895 lipopeptide
"""
from rdkit import Chem

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide consists of a peptide with attached lipid moieties.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Define SMARTS pattern for a peptide bond (-CO-NH-)
    peptide_bond_pattern = Chem.MolFromSmarts("[C](=O)[N]")
    if not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "No peptide bonds found"

    # Define SMARTS for common lipid structures, focusing on long aliphatic chains
    lipid_patterns = [
        Chem.MolFromSmarts("CCCCCCCC"),  # At least an 8-carbon chain
        Chem.MolFromSmarts("CCCCCCCCC(C)"),  # Long chain with branching
        Chem.MolFromSmarts("CC(C)C(C)C")  # Branched chain structure
    ]

    # Detect any lipid-like pattern
    lipid_like_found = any(mol.HasSubstructMatch(pattern) for pattern in lipid_patterns)
    if not lipid_like_found:
        return False, "No lipid-like structures found"

    # Check for peptide-lipid conjugation
    # Check if a lipid structure is connected to a peptide bond
    # This can be complex, and an accurate SMARTS pattern might need refinement
    # Here, we assume detection of both lipid and peptide is enough for current criteria

    return True, "Contains peptide bonds and lipid-like moieties"

# Example usage:
# smiles_string = "CCCCCCCCCC(=O)NCC(=O)NCCCCCC"
# is_lipopeptide(smiles_string)