"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: monoterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is any terpenoid derived from a monoterpene (C10 skeleton),
    including compounds in which the C10 skeleton has been rearranged or modified
    by the removal of one or more skeletal atoms (generally methyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Monoterpenoids are derived from monoterpenes (C10 skeleton), possibly modified
    if num_carbons < 7 or num_carbons > 10:
        return False, f"Number of carbon atoms ({num_carbons}) not consistent with monoterpenoids (should be between 7 and 10)"

    # Count the number of oxygen atoms
    num_oxygens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if num_oxygens == 0:
        return False, "No oxygen atoms present, unlikely to be a terpenoid"

    # Check for terpenoid functional groups
    functional_group_patterns = {
        "alcohol": "[OX2H]",                # Hydroxyl group
        "ketone": "[CX3](=O)[#6]",          # Ketone group
        "aldehyde": "[CX3H1](=O)[#6]",      # Aldehyde group
        "carboxylic_acid": "C(=O)[OH1]",    # Carboxylic acid group
        "ester": "C(=O)O[#6]",              # Ester group
        "ether": "[OD2]([#6])[#6]",         # Ether group
    }

    has_functional_group = False
    for fg_name, smarts in functional_group_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_functional_group = True
            break

    if not has_functional_group:
        return False, "No terpenoid functional groups detected"

    # Monoterpenoids may be cyclic or acyclic
    # Record if molecule is cyclic
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    # Check for branching (common in monoterpenoids)
    num_chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))

    # Compile reasons for classification
    reason = f"Contains {num_carbons} carbons, {num_oxygens} oxygens, "
    reason += f"{'cyclic' if num_rings > 0 else 'acyclic'} structure, "
    reason += f"{'with' if has_functional_group else 'without'} terpenoid functional groups"

    # If the molecule fits the criteria, classify as monoterpenoid
    return True, reason