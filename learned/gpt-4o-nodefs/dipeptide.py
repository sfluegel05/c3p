"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is likely to be a dipeptide based on its SMILES string.
    A dipeptide is generally composed of two amino acid units, possibly within more complex connections.
    
    It is vital to recognize varying structures, allowing some flexibility in peptide bond counts,
    provided there are recognizably two amino acid-like structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dipeptide variant, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for identifying peptide bonds (-CO-NH- linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Allow for one or two peptide bonds to account for variability in structure
    if len(peptide_bond_matches) < 1:
        return False, f"Expected at least 1 peptide bond, found {len(peptide_bond_matches)}"

    # SMARTS pattern for a more generalized amino acid core
    amino_acid_core_pattern = Chem.MolFromSmarts("[NX3][CX4]C(=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_core_pattern)
    
    # Check for at least two amino acid-like structures
    if len(amino_acid_matches) < 2:
        return False, f"Expected at least 2 amino acid residues, found {len(amino_acid_matches)}"

    return True, "Molecule likely contains two amino acid units that form a dipeptide structure"