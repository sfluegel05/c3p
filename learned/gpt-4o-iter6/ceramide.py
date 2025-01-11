"""
Classifies: CHEBI:17761 ceramide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    Ceramides are characterized by an amide linkage with fatty acids, a sphingoid base, and typically have hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ceramide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Look for amide linkage pattern
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide linkage found"

    # Check for long fatty acid chain, typically C14-C26
    carbon_chain_pattern = Chem.MolFromSmarts("[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]~[C]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found"

    # Check for sphingoid base, nitrogen must be linked in sphingoid region
    sphingoid_base_pattern = Chem.MolFromSmarts("NC[C@H](O)C")
    if not mol.HasSubstructMatch(sphingoid_base_pattern):
        return False, "No sphingoid base structure found"

    # Check for hydroxyl groups at positions consistent with ceramides
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "Missing hydroxyl group, which is common in ceramides"

    # Ensure overall structure supports ceramide classification
    n_atoms = mol.GetNumAtoms()
    n_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if n_nitrogens < 1 or n_carbons > n_atoms - n_oxygens - n_nitrogens:
        return False, "Structure inconsistent with ceramide"

    return True, "Structure contains features consistent with ceramide classification"