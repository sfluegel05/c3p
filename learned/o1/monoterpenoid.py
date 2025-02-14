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

    # Define SMARTS patterns for isoprene unit
    isoprene_pattern = Chem.MolFromSmarts('C=C(C)C')  # Isoprene unit

    # Find the number of isoprene units
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    num_isoprene_units = len(isoprene_matches)

    if num_isoprene_units < 2:
        # Monoterpenoids are derived from two isoprene units
        return False, f"Contains {num_isoprene_units} isoprene units, less than required for monoterpenoid"

    # Define SMARTS patterns for common monoterpene skeletons
    monoterpene_patterns = [
        "C1CCC(C=C1)C(C)C",        # Myrcene-like structures
        "C1=CC=CC=C1C(C)C",        # Limonene-like structures
        "C1=C(C)CCCC1C",           # Pinene-like structures
        "CC(C)C1CCC(C=C1)O",       # Menthol-like structures
        "CC(=C)C1CCC(C1)O",        # Pulegone-like structures
    ]

    has_monoterpene_skeleton = False
    for smarts in monoterpene_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_monoterpene_skeleton = True
            break

    if not has_monoterpene_skeleton:
        return False, "No monoterpene skeleton detected"

    # Check for terpenoid functional groups
    # Functional groups that might be present in monoterpenoids
    functional_group_patterns = {
        "alcohol": "[OX2H]",                # Hydroxyl group
        "ketone": "[CX3](=O)[#6]",          # Ketone group
        "aldehyde": "[CX3H1](=O)[#6]",      # Aldehyde group
        "ether": "[OD2]([#6])[#6]",         # Ether group
        "carboxylic_acid": "C(=O)[OH]",     # Carboxylic acid group
        "ester": "C(=O)O[#6]",              # Ester group
        "epoxide": "[C;R][O;R][C;R]",       # Epoxide ring
        "thiol": "[SX2H]",                  # Thiol group
        "double_bond": "C=C",               # Double bond
    }

    has_functional_group = False
    for fg_name, smarts in functional_group_patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            has_functional_group = True
            break

    if not has_functional_group:
        return False, "No terpenoid functional groups found"

    return True, "Contains monoterpene skeleton and terpenoid functional groups consistent with a monoterpenoid"