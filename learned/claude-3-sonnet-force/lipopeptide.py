"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: CHEBI:51359 lipopeptide
A lipopeptide is a compound consisting of a peptide with attached lipid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide backbone
    peptide_pattern = Chem.MolFromSmarts("[$([NX3H2,NX4H3+]),$([NX3H1,NX4H2+]([NX3H2,NX4H3+])([NX3H1,NX4H2+]))][CX3](=[OX1])[NX3H2,NX4H3+]")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    if not peptide_matches:
        return False, "No peptide backbone found"

    # Check for lipid chains
    lipid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    lipid_matches = mol.GetSubstructMatches(lipid_pattern)
    if not lipid_matches:
        return False, "No lipid chain found"

    # Check for ester or amide bonds between peptide and lipid
    bond_pattern = Chem.MolFromSmarts("[NX3H1,NX4H2+](-[CX3](=[OX1])-[OX2]-[CX4])")
    bond_matches = mol.GetSubstructMatches(bond_pattern)
    if not bond_matches:
        return False, "No ester or amide bond between peptide and lipid"

    # Check molecular weight - lipopeptides typically >500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for lipopeptide"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if c_count < 10 or o_count < 3 or n_count < 3:
        return False, "Insufficient atoms for lipopeptide"

    return True, "Contains both peptide backbone and lipid chain linked by ester or amide bond"