"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is typically defined as a peptide containing 2-20 amino acids.

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
    
    # Check for presence of peptide bonds (-C(=O)N-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    num_peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))
    
    if num_peptide_bonds < 1:
        return False, "No peptide bonds identified"
    if num_peptide_bonds >= 20:
        return False, f"Too many peptide bonds for an oligopeptide: {num_peptide_bonds}"
    
    # Assess if there are additional amino acids or modified residues
    amino_acid_patterns = [
        Chem.MolFromSmarts("[C@@H](N)C(=O)O"),  # Alpha amino acids
        Chem.MolFromSmarts("[C@H](N)C(=O)O"),   # Non-chiral centers can also appear in peptide sequences
        Chem.MolFromSmarts("N[C@@H](C)C(=O)O")  # Non-standard amino acids often derive from standard ones
    ]
    
    amino_acid_match = any(
        mol.HasSubstructMatch(pattern) for pattern in amino_acid_patterns
    )
    
    if not amino_acid_match:
        return False, "Missing recognizable amino acid residues"
    
    return True, "Oligopeptide structure identified with appropriate number of peptide bonds"