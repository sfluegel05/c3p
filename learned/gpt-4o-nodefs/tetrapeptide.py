"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide consists of exactly four amino acids linked by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a peptide bond SMILES pattern, focusing especially on 
    # N-C(=O) connections which are typical for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("N[C](=O)C")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Counts the peptide bonds and ensures there are exactly 3 (linking 4 residues)
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3"

    # Define a complete amino acid pattern matching an N-bound to an alpha-C and C=O
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4](C)[CX3](=O)[O,N]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # Each amino acid should have an N-C-C(=O)-Group pattern, counting 4 is crucial
    if len(amino_acid_matches) != 4:
        return False, f"Found {len(amino_acid_matches)} amino acid backbones, need exactly 4"

    return True, "Contains exactly four amino acids linked by three peptide bonds"