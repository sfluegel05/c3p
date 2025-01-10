"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is typically defined as a peptide containing a small number of amino acids, often between 2-20.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try to identify peptide bonds (amide linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    num_peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))
    
    if num_peptide_bonds < 1:
        return False, "No peptide bonds identified"
    
    # Check if the number of peptide bonds exceeds typical oligopeptide range
    if num_peptide_bonds > 25:
        return False, f"Excessive peptide bonds for an oligopeptide: {num_peptide_bonds}"
    
    # Recognize standard and modified amino acid residues
    amino_acid_patterns = [
        Chem.MolFromSmarts("[C@@H](N)C(=O)O"),  # L-alpha amino acids
        Chem.MolFromSmarts("[C@H](N)C(=O)O"),   # D-alpha amino acids
        Chem.MolFromSmarts("OC(=O)[C@@H](N)C"), # Alternative common structures
        Chem.MolFromSmarts("[N]C(=O)C")         # General amides for unusual side-chains
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in amino_acid_patterns):
        return False, "Missing recognizable amino acid residues"
    
    # Check for cyclic structures
    if mol.GetRingInfo().IsAtomInRingOfSize(3):
        return False, "Detected a very small ring, likely not an oligopeptide"
    
    # Additional criteria may include molecular weight checks, specific atom counts, etc.
    
    return True, "Oligopeptide structure identified with recognized peptide bonds and amino acid patterns"