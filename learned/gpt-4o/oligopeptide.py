"""
Classifies: CHEBI:25676 oligopeptide
"""
from rdkit import Chem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is defined as a sequence of amino acids connected by peptide bonds,
    typically with a relatively small number of amino acids (less than 20).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an oligopeptide, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string to create an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a peptide bond: NCC(=O)
    peptide_bond_pattern = Chem.MolFromSmarts("NCC(=O)")
    if peptide_bond_pattern is None:
        return False, "Unable to create SMARTS pattern for peptide bond"

    # Find matches for the peptide bond pattern
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # Check if the number of peptide bonds corresponds to an oligopeptide
    if num_peptide_bonds >= 1 and num_peptide_bonds < 20:
        return True, f"Contains {num_peptide_bonds} peptide bonds, consistent with an oligopeptide"
    elif num_peptide_bonds >= 20:
        return False, f"Contains {num_peptide_bonds} peptide bonds, too many for an oligopeptide"
    else:
        return False, "No peptide bonds found, not an oligopeptide"

# Example usage
smiles_example = "O=C(N[C@@H](CCC(O)=O)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC=1NC=NC1)CC=2NC=NC2"
is_oligopeptide(smiles_example)