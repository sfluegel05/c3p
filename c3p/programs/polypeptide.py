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

    # Define a broader range of amino acid residue patterns
    # Common peptide bond: -C(=O)-NH-
    peptide_bond_pattern = Chem.MolFromSmarts('C(=O)N')

    # Count the number of peptide bonds as a proxy for the number of amino acid residues
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # For polypeptides with 10+ residues, there should be at least 9 peptide bonds
    if num_peptide_bonds >= 9:
        return True, f"Contains {num_peptide_bonds + 1} amino acid residues, classifying as polypeptide"
    else:
        return False, f"Contains {num_peptide_bonds + 1} amino acid residues, not enough for a polypeptide"