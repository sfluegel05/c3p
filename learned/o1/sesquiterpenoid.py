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
    # Sesquiterpenoids typically have around 15 carbons, but can vary due to rearrangements or loss of methyl groups
    if c_count < 13 or c_count > 17:
        return False, f"Carbon count is {c_count}, which is not typical for a sesquiterpenoid (13-17 carbons including possible modifications)"

    # Check for presence of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structures found, which is uncommon for sesquiterpenoids"

    # Check for common terpenoid functional groups
    # (hydroxyls, ketones, aldehydes, carboxylic acids, esters, ethers, epoxides, double bonds)
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

    # Check that molecule is mostly composed of C, H, O, N, S atoms
    allowed_atomic_nums = {1, 6, 7, 8, 16}  # H, C, N, O, S
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains heteroatom {atom.GetSymbol()} not typical in sesquiterpenoids"

    # If all checks pass, classify as sesquiterpenoid
    return True, "Molecule meets criteria for a sesquiterpenoid (carbon count, functional groups, ring structures)"