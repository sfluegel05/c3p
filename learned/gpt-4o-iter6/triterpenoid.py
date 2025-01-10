"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is characterized by a core C30 structure derived from a triterpene,
    which can include rearrangements or modifications by removal of skeletal atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (25 <= c_count <= 35):  # Assuming some flexibility around 30 carbons
        return False, f"Carbon count {c_count} is not in the typical range for triterpenoids"

    # Check for polycyclic structure, common in terpenoids
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 4:
        return False, "Too few rings, typically polycyclic structures are seen in triterpenoids"

    # Check for isoprene units
    isoprene_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene-like units found"

    # Check for common functional groups like hydroxyl, carbonyl, etc.
    # These checks ensure that we have functional group diversity typical of triterpenoids
    if not any(Chem.MolFromSmarts(sm).HasSubstructMatch(mol) for sm in ["C=O", "O"]):
        return False, "No characteristic functional groups (e.g., carbonyl or hydroxyl) found"

    return True, "Structure conforms with general criteria for triterpenoids"