"""
Classifies: CHEBI:48030 tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide contains four amino-acid residues connected by peptide linkages.

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

    # Define revised SMARTS pattern for peptide bonds, considering more diverse connectivity
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[C]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3 for tetrapeptide"

    # Check for the presence of four amino acid residues by expanding scope
    aa_residue_pattern = Chem.MolFromSmarts("[N](C[CX4])[CX3](=O)[C;R0]")
    aa_residue_matches = mol.GetSubstructMatches(aa_residue_pattern)
    
    if len(aa_residue_matches) != 4:
        return False, f"Detected {len(aa_residue_matches)} amino acid residues, need exactly 4"
    
    return True, "Contains four amino-acid residues connected by three peptide bonds, valid tetrapeptide structure"