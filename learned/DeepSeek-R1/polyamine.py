"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:53371 polyamine
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmarts

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine contains two or more amino groups (-NH2, -NH-, etc.) in an organic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Predefined SMARTS patterns for exclusion groups
    nitro_pattern = MolFromSmarts('[NX3](=[OX1])=[OX1]')  # Nitro groups
    nitrile_pattern = MolFromSmarts('[C]#[N]')            # Nitriles
    sulfonamide_pattern = MolFromSmarts('[SX4](=[OX1])(=[OX1])([#6])[NX3]')  # Sulfonamides
    amide_pattern = MolFromSmarts('[CX3](=[OX1])[NX3H0,H1]')  # Amides excluding those with O substituents

    amine_count = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Only process nitrogen atoms

        # Check if nitrogen is part of any exclusion group
        is_excluded = False
        for pattern in [nitro_pattern, nitrile_pattern, sulfonamide_pattern, amide_pattern]:
            if mol.HasSubstructMatch(pattern) and atom.HasSubstructMatch(pattern):
                is_excluded = True
                break
        if is_excluded:
            continue

        # Check if nitrogen is bonded to at least one carbon (organic compound requirement)
        has_carbon = any(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
        if has_carbon:
            amine_count += 1

    if amine_count >= 2:
        return True, f"Contains {amine_count} amine groups"
    return False, f"Only {amine_count} amine groups found (requires ≥2)"