"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide consists of three amino acid residues connected by peptide linkages.

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
    
    # Define the peptide bond pattern: N-C(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Ensure exactly 2 peptide bonds since tripeptides consist of three residues connected by two peptide bonds
    if len(peptide_bond_matches) != 2:
        return False, f"Expected 2 peptide bonds, found {len(peptide_bond_matches)}"
    
    # Check for presence of 3 amino acid residues
    # Amino acids typically have a structure like: N-C-C(=O)
    amino_acid_pattern = Chem.MolFromSmarts("N-C-C(=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Count unique alpha carbon atoms bonded to amino groups and carboxyl groups
    n_amino_acids = len(set(match[1] for match in amino_acid_matches))  # Index 1 is the central C in N-C-C(=O)

    # Validate exactly 3 amino acids
    if n_amino_acids != 3:
        return False, f"Molecule spans {n_amino_acids} amino acids, should have 3"
    
    return True, "Molecule is a tripeptide with three amino acids connected by peptide linkages"