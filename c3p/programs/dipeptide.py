"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is a molecule that contains two amino acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an alpha carbon in a peptide chain
    # Alpha carbon connected to nitrogen, carbonyl carbon, and side chain (R group)
    alpha_carbon_smarts = '[CX4H1]([#6])[CX3](=O)[#7]'
    alpha_carbon = Chem.MolFromSmarts(alpha_carbon_smarts)
    if alpha_carbon is None:
        return False, "Invalid alpha carbon SMARTS pattern"

    # Find all alpha carbons
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon)
    num_alpha_carbons = len(alpha_carbons)

    # Define the SMARTS pattern for peptide bond (amide bond between carbonyl carbon and nitrogen)
    peptide_bond_smarts = 'C(=O)N'
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond is None:
        return False, "Invalid peptide bond SMARTS pattern"

    # Find all peptide bonds
    peptide_bonds = mol.GetSubstructMatches(peptide_bond)
    num_peptide_bonds = len(peptide_bonds)

    # Check if the molecule has exactly two amino acid residues and one peptide bond
    if num_alpha_carbons == 2 and num_peptide_bonds == 1:
        return True, "Contains two amino acid residues connected by a peptide bond"
    else:
        return False, f"Found {num_alpha_carbons} amino acid residues and {num_peptide_bonds} peptide bonds"