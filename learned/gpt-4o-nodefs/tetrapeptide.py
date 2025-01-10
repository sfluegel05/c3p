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

    # Define a more selective peptide bond SMARTS pattern
    # This pattern now includes the nitrogen followed by the carbonyl-carbon that is bonded directly to the C-alpha of the amino acid
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Counts the peptide bonds and ensures there are exactly 3
    if len(peptide_bond_matches) != 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need exactly 3"

    # Further checks: Verify there are exactly four alpha-amino acid backbones
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])")
    amino_acid_matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)

    # Each alpha amino acid should contain this pattern, count to see if there are 4
    if len(amino_acid_matches) != 4:
        return False, f"Found {len(amino_acid_matches)} alpha-amino groups, need exactly 4"

    return True, "Contains exactly four amino acids linked by three peptide bonds"