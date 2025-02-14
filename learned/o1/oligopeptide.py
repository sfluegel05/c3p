"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: oligopeptide

An oligopeptide is a peptide containing a relatively small number of amino acids.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide based on its SMILES string.
    An oligopeptide is defined as a peptide consisting of 2 to 20 amino acid residues linked by peptide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oligopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify peptide (amide) bonds between amino acids
    # Peptide bond SMARTS pattern: N-C(=O)
    peptide_bond_smarts = "[NX3][CX3](=O)[#6]"
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond)
    
    # Count peptide bonds
    num_peptide_bonds = len(peptide_bond_matches)
    if num_peptide_bonds == 0:
        return False, "No peptide bonds found"

    # Estimate number of amino acid residues (number of peptide bonds + 1)
    num_residues = num_peptide_bonds + 1
    
    # Check if the number of residues falls within the oligopeptide range
    if 2 <= num_residues <= 20:
        return True, f"Molecule has {num_residues} amino acid residues, classified as an oligopeptide"
    else:
        return False, f"Molecule has {num_residues} amino acid residues, not within oligopeptide range (2-20)"