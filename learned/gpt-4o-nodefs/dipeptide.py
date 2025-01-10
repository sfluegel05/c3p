"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide consists of two amino acid residues connected by a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a peptide bond
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    
    if len(peptide_bond_matches) != 1:
        return False, f"Expected 1 peptide bond, found {len(peptide_bond_matches)}"

    # SMARTS pattern for amino acid residues (basic pattern)
    amino_acid_backbone_pattern = Chem.MolFromSmarts("[NX3][CX4][CX3](=[OX1])[#6]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_backbone_pattern)

    if len(amino_acid_matches) != 2:
        return False, f"Expected 2 amino acid residues, found {len(amino_acid_matches)}"

    return True, "Molecule contains exactly 2 amino acid residues connected by a peptide bond, forming a dipeptide"