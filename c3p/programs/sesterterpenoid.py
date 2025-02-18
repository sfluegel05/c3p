"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from a sesterterpene, often featuring rearranged or modified C25 backbones.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Isoprene units are 5 carbons; sesterterpenoids typically have 25-carbon backbones
    # Check if the carbon count is in a reasonable range for a sesterterpenoid (near 25)
    if c_count < 20 or c_count > 35:
        return False, f"Carbon count of {c_count} is outside typical range for sesterterpenoids"

    # Check for functional groups common in sesterterpenoids: hydroxyl (OH), ester, and ketone
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    ketone_pattern = Chem.MolFromSmarts("C(=O)[CX4]")
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CX4]")

    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    has_ketone = mol.HasSubstructMatch(ketone_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)

    # Typically, sesterterpenoids will have these functional groups due to their biosynthetic derivations
    if not (has_hydroxyl or has_ketone or has_ester):
        return False, "Missing common functional groups: hydroxyl, ketone, or ester"

    # Presence of isoprene components - although rearranged, some patterns should be similar
    has_terpenoid_signature = False
    isoprene_pattern = Chem.MolFromSmarts("[CX4]1[CX2]=[CX1][CX2][CX4]1")
    if mol.HasSubstructMatch(isoprene_pattern):
        has_terpenoid_signature = True

    if not has_terpenoid_signature:
        return False, "Does not have clear terpenoid signature"

    return True, "Contains characteristics typical of a sesterterpenoid including carbon count and functional groups"