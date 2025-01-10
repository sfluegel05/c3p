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
        bool: True if the molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for peptide bond pattern (C(=O)N)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Count peptide bonds
    n_peptide_bonds = len(peptide_bond_matches)
    
    # Polypeptides must have at least 10 peptide bonds (for ten or more residues)
    if n_peptide_bonds >= 10:
        return True, f"Contains {n_peptide_bonds} peptide bonds, indicating a polypeptide"
    else:
        return False, f"Contains {n_peptide_bonds} peptide bonds, fewer than required for a polypeptide"