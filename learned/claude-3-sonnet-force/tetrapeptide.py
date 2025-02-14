"""
Classifies: CHEBI:48030 tetrapeptide
"""
"""
Classifies: CHEBI:25104 tetrapeptide
A tetrapeptide is any molecule that contains four amino-acid residues connected by peptide linkages.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    A tetrapeptide is any molecule that contains four amino-acid residues connected by peptide linkages.

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

    # Check for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[NX3]")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 3:
        return False, f"Found {len(peptide_bond_matches)} peptide bonds, need at least 3 for tetrapeptide"

    # Check for amino acid residues
    amino_acid_pattern = Chem.MolFromSmarts("[NX3][CX4H]([CH2])[CX3](=[OX1])[OX2H,OX1-,N]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if len(set(amino_acid_matches)) != 4:
        return False, f"Found {len(set(amino_acid_matches))} amino acid residues, need exactly 4 for tetrapeptide"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 800:
        return False, "Molecular weight outside typical range for tetrapeptides"

    return True, "Molecule contains four amino acid residues connected by peptide linkages"