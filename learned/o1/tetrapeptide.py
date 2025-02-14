"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: tetrapeptide
"""
from rdkit import Chem

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is a molecule containing four amino acid residues connected by peptide linkages.

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

    # Find peptide bonds (amide bonds between amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)

    num_peptide_bonds = len(peptide_bond_matches)

    if num_peptide_bonds != 3:
        return False, f"Found {num_peptide_bonds} peptide bonds, expected 3 for tetrapeptide"

    # Identify alpha carbons (Cα) - carbon attached to both an amine and a carboxyl group
    alpha_carbon_pattern = Chem.MolFromSmarts("[C;!R][C](=O)[N]")  # Non-ring carbon connected to C(=O)N
    alpha_carbons = mol.GetSubstructMatches(alpha_carbon_pattern)
    num_alpha_carbons = len(alpha_carbons)

    if num_alpha_carbons != 4:
        return False, f"Found {num_alpha_carbons} alpha carbons, expected 4 for tetrapeptide"

    # Verify that the peptide bonds connect the amino acids in sequence
    # Check that the molecule is mostly linear (no significant branching)
    n_chains = Chem.GetMolFrags(mol, asMols=False, sanitizeFrags=False)
    if len(n_chains) > 1:
        return False, "Molecule is fragmented, expected a single continuous peptide chain"

    # Optionally, check for terminal groups (N-terminus and C-terminus)
    # N-terminus amine group
    n_term_pattern = Chem.MolFromSmarts("N[^C(=O)]")  # Nitrogen not attached to C=O
    n_term_matches = mol.GetSubstructMatches(n_term_pattern)
    if len(n_term_matches) == 0:
        return False, "No N-terminal amine group found"

    # C-terminus carboxyl group
    c_term_pattern = Chem.MolFromSmarts("C(=O)[O;H1,H0-]")  # C=O attached to OH or O-
    c_term_matches = mol.GetSubstructMatches(c_term_pattern)
    if len(c_term_matches) == 0:
        return False, "No C-terminal carboxyl group found"

    # Ensure there are no additional peptide bonds or amino acids
    # Total expected atoms in tetrapeptide backbone: (4 residues) * (N-Cα-C') = 12 atoms
    # plus side chains and terminal groups
    # For simplicity, assuming the molecule doesn't contain additional peptide bonds
    return True, "Molecule is a tetrapeptide with 4 amino acid residues connected via peptide bonds"