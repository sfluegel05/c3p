"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from a sesterterpene, showcasing a modified 25-carbon backbone, rearranged 
    or involving ring structures and diverse functional groups.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a sesterterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Sesterterpenoids often have ~25 carbons; allow for minor modifications
    if c_count < 20 or c_count > 30:
        return False, f"Carbon count of {c_count} is outside the range typical for sesterterpenoids"

    # Check for rings; sesterterpenoids frequently have multiple rings
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count < 2:
        return False, "Contains less than 2 ring structures typical for sesterterpenoids"

    # Look for common terpenoid structure indicator: isoprene (C5) like units
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "Lacks recognized terpenoid C5-like structural units"

    # Functional group diversity; check for heteroatom presence (N, O, S)
    heteroatoms_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8, 16])
    if heteroatoms_count < 1:
        return False, "Insufficient heteroatoms indicating lack of functional group diversity common in sesterterpenoids"

    # If all checks are satisfied, classify as sesterterpenoid
    return True, "Molecule has characteristics typical of a sesterterpenoid: modified backbone, functional diversity, and sufficient ring structures"