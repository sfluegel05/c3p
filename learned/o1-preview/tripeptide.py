"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is an oligopeptide that consists of three amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tripeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for peptide bond (amide bond in backbone)
    peptide_bond_smarts = "[C;$(C(=O))][N;H1,H2]"
    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)

    # Find peptide bonds
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)

    # Tripeptides have two peptide bonds connecting three amino acid residues
    if num_peptide_bonds != 2:
        return False, f"Found {num_peptide_bonds} peptide bonds, expected 2 for tripeptide"

    # Identify amino acid residues by finding alpha-carbon patterns
    alpha_carbon_smarts = "[N;H1,H2][C@@H](*)C(=O)"
    alpha_carbon_pattern = Chem.MolFromSmarts(alpha_carbon_smarts)
    residues = mol.GetSubstructMatches(alpha_carbon_pattern)
    num_residues = len(residues)

    # Allow for N-terminal modifications by checking if residues >= 3
    if num_residues < 3:
        return False, f"Found {num_residues} amino acid residues, expected at least 3 for tripeptide"

    # Additional check: ensure that residues are connected via peptide bonds
    # Build a graph of residues connected by peptide bonds
    residue_atoms = set()
    for match in residues:
        residue_atoms.update(match)

    # Check connectivity via peptide bonds
    peptide_bond_atoms = []
    for match in peptide_bonds:
        peptide_bond_atoms.extend(match)

    connected_residues = residue_atoms.intersection(peptide_bond_atoms)

    if len(connected_residues) < 2:
        return False, "Amino acid residues are not connected via peptide bonds"

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"