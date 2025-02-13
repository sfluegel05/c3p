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

    # Pattern specifying C(α)-N in -CONH- linkage, part of peptide backbone
    alpha_carbon_pattern = Chem.MolFromSmarts("N[C@H]")
    # Peptide bond pattern: -C(=O)N-
    peptide_pattern = Chem.MolFromSmarts("C(=O)N")
    
    # Check for valid patterns
    if alpha_carbon_pattern is None or peptide_pattern is None:
        return False, "Pattern creation failed"

    # Find matching functionalities, crucial to amino acids
    alpha_matches = mol.GetSubstructMatches(alpha_carbon_pattern)
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)

    # Check if it has three alpha carbon matches indicating three residues
    if len(alpha_matches) == 3 and len(peptide_matches) >= 2:
        return True, "Contains three amino acid residues connected by peptide linkages"
    else:
        return False, f"Contains {len(alpha_matches)} alpha carbons and {len(peptide_matches)} peptide bonds, insufficient or excess for a tripeptide"