"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Lipinski

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is defined as a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least 10 amino acid residues
    amino_acid_pattern = Chem.MolFromSmarts('[NH2][CX4H]([CH3])([CH3])[CX4H]([NH2])[CX3](=[OH])[OH]')
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(amino_acid_matches) < 10:
        return False, f"Found {len(amino_acid_matches)} amino acid residues, but at least 10 are required for a polypeptide"

    # Check for linear or cyclic connectivity of amino acid residues
    linear_pattern = Chem.MolFromSmarts('[NH2][CX4H]([CH3])([CH3])[CX4H]([NH2])[CX3](=[OH])[NH][CX4H]([CH3])([CH3])[CX4H]([NH2])[CX3](=[OH])')
    cyclic_pattern = Chem.MolFromSmarts('[NH2][CX4H]([CH3])([CH3])[CX4H]([NH2])[CX3](=[OH])[NH]1[CX4H]([CH3])([CH3])[CX4H]([NH2])[CX3](=[OH])[NH][CX4H]([CH3])([CH3])[CX4H]([NH2])[CX3](=[OH])[NH]1')
    if not mol.HasSubstructMatch(linear_pattern) and not mol.HasSubstructMatch(cyclic_pattern):
        return False, "Amino acid residues are not connected in a linear or cyclic fashion"

    # Check for amide bonds and peptide bonds
    amide_bond_pattern = Chem.MolFromSmarts('[NH][CX3](=[OH])[CX3]=[OH]')
    peptide_bond_pattern = Chem.MolFromSmarts('[NH][CX3](=[OH])[NH][CX4H]([CH3])([CH3])[CX4H]([NH2])[CX3](=[OH])')
    if not mol.HasSubstructMatch(amide_bond_pattern) or not mol.HasSubstructMatch(peptide_bond_pattern):
        return False, "Missing characteristic amide or peptide bonds"

    # Check for conformance with Lipinski's rule of 5
    if not Lipinski.CompliesWithRuleOfFive(mol):
        return False, "Does not comply with Lipinski's rule of 5, which is typically expected for polypeptides"

    return True, "Contains at least 10 amino acid residues connected in a linear or cyclic fashion, with characteristic amide and peptide bonds"