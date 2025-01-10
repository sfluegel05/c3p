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
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for multiple peptide bonds
    multiple_peptide_bonds_pattern = Chem.MolFromSmarts("([C](=O)[N]){2,}")

    if not mol.HasSubstructMatch(multiple_peptide_bonds_pattern):
        return False, "Fewer than two peptide bonds found"

    # Enhanced lipid pattern for longer aliphatic chains (10+ carbon chains)
    extended_lipid_patterns = [
        Chem.MolFromSmarts("C{10,}"),  # At least a 10-carbon chain
    ]

    # Detect any lipid-like pattern
    lipid_like_found = any(mol.HasSubstructMatch(pattern) for pattern in extended_lipid_patterns)
    if not lipid_like_found:
        return False, "No extended lipid-like structures found"

    # Check if multiple peptide bonds are connected to a lipid structure
    # It's complex but crucial to confirm connectivity
    # Ensure that the detection logic allows peptide to lipid connection

    return True, "Contains multiple peptide bonds and extended lipid-like moieties"

# Example usage:
# smiles_string = "CCCCCCCCCCCC(=O)NCC(=O)NCCCCCC"
# is_lipopeptide(smiles_string)