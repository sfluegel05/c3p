"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid is a terpenoid derived from a sesquiterpene, built from three isoprene units (C5 units),
    typically containing around 15 carbons, but may have rearrangements or modifications by the removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # Count other heteroatoms (excluding hydrogen)
    heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1,6])

    # Sesquiterpenoids typically have around 15 carbons
    if c_count < 12 or c_count > 15:
        return False, f"Carbon count is {c_count}, which is not typical for a sesquiterpenoid (12-15 carbons including possible modifications)"

    # Check for terpenoid functional groups (hydroxyls, ketones, aldehydes, carboxyls, esters, epoxides)
    functional_groups = {
        'hydroxyl': '[OX2H]',                         # Alcohol
        'ketone': '[CX3](=O)[#6]',                   # Ketone
        'aldehyde': '[CX3H1](=O)',                   # Aldehyde
        'carboxylic acid': '[CX3](=O)[OX2H1]',       # Carboxylic acid
        'ester': '[CX3](=O)[OX2H0][#6]',             # Ester
        'ether': '[OD2]([#6])[#6]',                  # Ether
        'epoxide': 'C1OC1',                          # Epoxide ring
        'double bond': 'C=C',                        # Double bond
    }
    fg_found = False
    for fg_name, fg_smarts in functional_groups.items():
        fg_pattern = Chem.MolFromSmarts(fg_smarts)
        if fg_pattern and mol.HasSubstructMatch(fg_pattern):
            fg_found = True
            break

    if not fg_found:
        return False, "No typical terpenoid functional groups found"

    # Check for isoprene units (C5 units)
    isoprene_pattern = Chem.MolFromSmarts('C(=C)C=C(C)C')
    if isoprene_pattern:
        isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
        if len(isoprene_matches) < 2:
            return False, f"Found {len(isoprene_matches)} isoprene units, less than expected for a sesquiterpenoid"
    else:
        return False, "Invalid isoprene SMARTS pattern"

    # Optional: Check for ring systems common in sesquiterpenoids (e.g., 5 or 6 membered rings)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structures found, which is uncommon for sesquiterpenoids"

    # If all checks pass, classify as sesquiterpenoid
    return True, "Molecule meets criteria for a sesquiterpenoid (carbon count, functional groups, isoprene units, ring structures)"