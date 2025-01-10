"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem

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

    # Define peptide bond pattern: amide bond between carbonyl carbon and nitrogen
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bond_matches)

    # For tripeptide, there should be exactly 2 peptide bonds
    if num_peptide_bonds != 2:
        return False, f"Found {num_peptide_bonds} peptide bonds, expected 2 for tripeptide"

    # Define amino acid residue pattern: backbone N-C(alpha)-C(=O)
    amino_acid_pattern = Chem.MolFromSmarts("N[C@@H]C(=O)")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    num_residues = len(amino_acid_matches)

    # Some amino acids might not match the pattern due to stereochemistry; try a more general pattern
    if num_residues == 0:
        amino_acid_pattern = Chem.MolFromSmarts("NCC(=O)")
        amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
        num_residues = len(amino_acid_matches)

    # For tripeptide, there should be exactly 3 amino acid residues
    if num_residues != 3:
        return False, f"Found {num_residues} amino acid residues, expected 3 for tripeptide"

    # Optionally, verify that peptide bonds connect the amino acids
    # For each peptide bond, check if it connects two amino acid residues
    connected_residues = 0
    for match in peptide_bond_matches:
        c_index = match[0]  # Carbonyl carbon
        n_index = match[1]  # Amide nitrogen
        c_atom = mol.GetAtomWithIdx(c_index)
        n_atom = mol.GetAtomWithIdx(n_index)

        # Check if the C and N atoms are part of amino acid backbones
        c_in_residue = False
        n_in_residue = False
        for residue_match in amino_acid_matches:
            if c_index in residue_match:
                c_in_residue = True
            if n_index in residue_match:
                n_in_residue = True

        if c_in_residue and n_in_residue:
            connected_residues += 1

    if connected_residues != 2:
        return False, "Peptide bonds do not properly connect amino acid residues"

    return True, "Molecule is a tripeptide with three amino acid residues connected by peptide bonds"