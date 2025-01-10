"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is composed of two amino acid residues joined by a single peptide bond.
    
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

    # SMARTS pattern for identifying a peptide bond (-CO-NH- linkage)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    # Check for exactly one peptide bond
    if len(peptide_bond_matches) != 1:
        return False, f"Expected 1 peptide bond, found {len(peptide_bond_matches)}"

    # SMARTS pattern for identifying common amino acid structures
    # A more relaxed pattern considering various forms and possible protective groups
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$([NX3][O,S])]C(=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)

    # Check for exactly two amino acid residues
    if len(amino_acid_matches) != 2:
        return False, f"Expected 2 amino acid residues, found {len(amino_acid_matches)}"

    # Additional verification could involve checking the side chains but for now, 
    # we check the structure that allows for variability.

    return True, "Molecule contains exactly 2 amino acid residues connected by a peptide bond, forming a dipeptide"