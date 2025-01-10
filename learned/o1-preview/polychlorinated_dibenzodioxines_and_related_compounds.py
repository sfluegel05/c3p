"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds
"""

from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl,
    or a related compound based on its SMILES string.

    The compound must contain one of the core structures (dibenzodioxin, dibenzofuran, or biphenyl)
    and must be polychlorinated or polybrominated (containing at least 3 halogen atoms).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polychlorinated dibenzodioxin or related compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structures using SMARTS patterns
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')  # Biphenyl core
    dibenzodioxin_pattern = Chem.MolFromSmarts('O1c2ccccc2Oc3ccccc13')  # Dibenzodioxin core
    dibenzofuran_pattern = Chem.MolFromSmarts('c1cc2oc3ccccc3cc2cc1')  # Dibenzofuran core

    # Check for core structures
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)
    has_dibenzodioxin = mol.HasSubstructMatch(dibenzodioxin_pattern)
    has_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)

    if not (has_biphenyl or has_dibenzodioxin or has_dibenzofuran):
        return False, "Does not contain dibenzodioxin, dibenzofuran, or biphenyl core structure"

    # Count Cl and Br atoms
    num_Cl = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl')
    num_Br = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Br')

    num_halogen = num_Cl + num_Br

    if num_halogen < 3:
        return False, f"Contains {num_halogen} halogen atoms; needs at least 3 to be polyhalogenated"

    # Check if halogens are attached to the aromatic rings
    halogen_on_ring = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in ['Cl', 'Br']:
            for neighbor in atom.GetNeighbors():
                if neighbor.IsInRing() and neighbor.GetIsAromatic():
                    halogen_on_ring = True
                    break
        if halogen_on_ring:
            break

    if not halogen_on_ring:
        return False, "Halogens are not attached to aromatic rings"

    return True, "Contains polyhalogenated dibenzodioxin, dibenzofuran, or biphenyl core"