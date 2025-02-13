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

    # Improved alpha carbon chiral centers with side chains
    # This pattern aims to be more flexible regarding alpha carbon configuration
    alpha_carbon_patterns = [
        Chem.MolFromSmarts("NC[C@H]"),
        Chem.MolFromSmarts("NC[C@@H]"),
        Chem.MolFromSmarts("N[C@]C"),
        Chem.MolFromSmarts("N[C@@]C")
    ]
    
    # Peptide bond pattern: -C(=O)N- (amide linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")

    # Check for valid patterns
    if peptide_bond_pattern is None or any(pat is None for pat in alpha_carbon_patterns):
        return False, "Pattern creation failed"

    # Finding matches
    alpha_matches = sum([len(mol.GetSubstructMatches(pat)) for pat in alpha_carbon_patterns])
    peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))

    # Tripeptides have 3 residues and exactly 2 peptide bonds
    if alpha_matches == 3 and peptide_bonds == 2:
        return True, "Contains three alpha carbon chiral centers with two peptide linkages indicating a tripeptide"
    else:
        return False, f"Contains {alpha_matches} alpha carbons and {peptide_bonds} peptide bonds, insufficient or excess for a tripeptide"