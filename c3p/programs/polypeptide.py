"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General amino acid residue pattern, including common side chains
    amino_acid_pattern = Chem.MolFromSmarts('N[C@@H](C(=O))')

    # Find all matches of the amino acid pattern
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # Count number of amino acids based on the pattern
    num_amino_acids = len(amino_acid_matches)

    if num_amino_acids >= 10:
        return True, f"Contains {num_amino_acids} amino acid residues, classifying as polypeptide"
    else:
        return False, f"Contains {num_amino_acids} amino acid residues, not enough for a polypeptide"