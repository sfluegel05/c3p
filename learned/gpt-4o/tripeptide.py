"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide that consists of three amino acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Peptide bond pattern: (-C(=O)N-)
    peptide_pattern = Chem.MolFromSmarts("C(=O)N")
    if peptide_pattern is None:
        return False, "Failed to create peptide pattern"

    # Find all peptide bonds in the molecule
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    
    # Check if there are exactly two peptide bonds
    if len(peptide_matches) == 2:
        return True, "Contains exactly two peptide bonds indicating a tripeptide"
    else:
        return False, f"Contains {len(peptide_matches)} peptide bonds, which does not indicate a tripeptide"