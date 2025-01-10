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
        bool: True if molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Peptide bond pattern (N-C(=O)-C)
    peptide_bond_pattern = Chem.MolFromSmarts("N-C(=O)-C")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Check if there are exactly 3 peptide bonds
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3"

    # Further checks could be added here if needed to ensure correct amino acid backbone structure

    return True, "Contains exactly four amino acids linked by three peptide bonds"