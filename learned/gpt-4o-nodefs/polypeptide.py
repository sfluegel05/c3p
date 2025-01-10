"""
Classifies: CHEBI:15841 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    Polypeptides are long chains of amino acids with peptide bonds.

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

    # Look for peptide bond pattern (-CO-NH-)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bond_matches) < 2:
        return False, f"Detected {len(peptide_bond_matches)} peptide bonds, need at least 2 for a polypeptide"

    # Check for characteristic long backbone sequences
    backbone_pattern = Chem.MolFromSmarts("N[C,C](C(=O))[N]")
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if len(backbone_matches) < 1:
        return False, "No characteristic polypeptide backbone found"

    # Check for sufficient molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for polypeptide"

    return True, "Contains peptide bonds and polypeptide backbone structure"