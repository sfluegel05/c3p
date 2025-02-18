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

    # Peptide bond SMARTS pattern: represents an amide bond typically found between amino acids
    peptide_bond_pattern = Chem.MolFromSmarts('C(=O)N')

    # Find all matches of the peptide bond pattern
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Estimate number of amino acids: because each peptide bond connects two amino acids
    num_amino_acids = len(peptide_bond_matches) + 1

    if num_amino_acids >= 10:
        return True, f"Contains {num_amino_acids} amino acid residues, classifying as polypeptide"
    else:
        return False, f"Contains {num_amino_acids} amino acid residues, not enough for a polypeptide"