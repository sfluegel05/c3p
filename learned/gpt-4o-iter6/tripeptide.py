"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as three amino-acid residues connected by two peptide linkages (amide bonds).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Correct pattern: three amide bonds (N-C(=O) linkage)
    amide_bond_pattern = Chem.MolFromSmarts("NC(=O)C")
    amide_matches = mol.GetSubstructMatches(amide_bond_pattern)

    if len(amide_matches) != 2:
        return False, f"Found {len(amide_matches)} amide bonds, need exactly 2 for tripeptide linkage"

    # Check for the presence of three alpha carbon centers (chiral centers next to amide)
    alpha_carbon_pattern = Chem.MolFromSmarts("[C@H]([N])")
    alpha_carbon_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    
    if not (2 <= len(alpha_carbon_matches) <= 3):  # Flexible to small variations in chirality
        return False, f"Expected 2 or 3 chiral centers, got {len(alpha_carbon_matches)}"
    
    # Additional checks for cyclicity if needed, etc., could be here

    return True, "Contains three amino-acid residues connected by two peptide linkages"