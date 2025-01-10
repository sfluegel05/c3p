"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

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
    
    # Improved peptide bond pattern considering branching and stereochemistry
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Count peptide bonds
    if len(peptide_bond_matches) < 2:
        return False, "Insufficient peptide bonds, need at least 2 to form a tripeptide"
    
    # Amino acid backbone pattern more inclusive
    amino_acid_backbone_pattern = Chem.MolFromSmarts("[N][CX3](=O)[CX4]")
    backbone_matches = mol.GetSubstructMatches(amino_acid_backbone_pattern)

    # Count amino acid residues (each matching the pattern once forms part of a single amino acid)
    n_amino_acids = len(backbone_matches)
    if n_amino_acids != 3:
        return False, f"Molecule has {n_amino_acids} amino acids, should have 3"
    
    return True, "Molecule is a tripeptide with three amino acids connected by peptide linkages"