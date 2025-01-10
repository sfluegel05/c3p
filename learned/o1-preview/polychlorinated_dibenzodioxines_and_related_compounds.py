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
    and have at least one chlorine or bromine atom attached directly to the aromatic carbons of the core.

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
    dibenzofuran_pattern = Chem.MolFromSmarts('c1cc2coc3ccccc3cc2cc1')  # Dibenzofuran core

    # Check for core structures
    core_patterns = []
    if mol.HasSubstructMatch(biphenyl_pattern):
        core_patterns.append(('biphenyl', biphenyl_pattern))
    if mol.HasSubstructMatch(dibenzodioxin_pattern):
        core_patterns.append(('dibenzodioxin', dibenzodioxin_pattern))
    if mol.HasSubstructMatch(dibenzofuran_pattern):
        core_patterns.append(('dibenzofuran', dibenzofuran_pattern))

    if not core_patterns:
        return False, "Does not contain dibenzodioxin, dibenzofuran, or biphenyl core structure"

    # Define halogen atoms
    halogens = ['Cl', 'Br']

    # Check for halogens attached to the core structures
    halogen_on_core = False
    for core_name, core_pattern in core_patterns:
        core_matches = mol.GetSubstructMatches(core_pattern)
        for match in core_matches:
            core_atom_indices = set(match)
            for atom_idx in core_atom_indices:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetSymbol() == 'C' and atom.GetIsAromatic():
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetIdx() not in core_atom_indices:
                            if neighbor.GetSymbol() in halogens:
                                halogen_on_core = True
                                break
                    if halogen_on_core:
                        break
            if halogen_on_core:
                break
        if halogen_on_core:
            break

    if not halogen_on_core:
        return False, "No halogens attached to core structure"

    return True, "Contains polychlorinated/polybrominated dibenzodioxin, dibenzofuran, or biphenyl core"