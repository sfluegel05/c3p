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

    # Define SMARTS pattern for peptide bonds, considering possible variations
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    if len(peptide_bond_matches) < 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3 for tetrapeptide"

    # Define SMARTS pattern for typical amino acid residues
    # This includes patterns for amino and carboxyl groups connected by alpha carbon
    # Cannot rely entirely on SMARTS here due to the variability of peptide side chains
    amide_linked = Chem.MolFromSmarts("[NX3][CX3](=O)[C;R0]")

    # Check for typical alpha carbon, amide/link patterns should be >= 4 chunks matching typical amino acid sections
    aa_residue_count = 0
    for residue in mol.GetSubstructMatches(amide_linked):
        aa_residue_count += 1
    
    # Validate that there are four amino acid residues
    if aa_residue_count != 4:
        return False, f"Detected {aa_residue_count} amino acid residues, need exactly 4"

    return True, "Contains four amino-acid residues connected by three peptide bonds, valid tetrapeptide structure"