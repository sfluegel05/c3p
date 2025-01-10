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
        bool: True if the molecule is a tetrapeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General peptide bond pattern: N-C(=O) linked to another N
    peptide_bond_pattern = Chem.MolFromSmarts("N-[C](=O)-N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    # Improve detection by checking continuity of peptide bonds (linked sequence)
    if len(peptide_bond_matches) < 3:
        return False, f"Contains {len(peptide_bond_matches)} peptide bonds, expected at least 3"

    # Define an improved pattern for amino acid residues in a peptide chain.
    amino_acid_pattern = Chem.MolFromSmarts("N[C](C)C(=O)")  # more specific to recognize standard linkage
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Check if we have exactly four residues
    if len(amino_acid_matches) != 4:
        return False, f"Detected {len(amino_acid_matches)} amino acid residues, expected exactly 4"

    return True, "Contains four amino-acid residues connected by peptide linkages"