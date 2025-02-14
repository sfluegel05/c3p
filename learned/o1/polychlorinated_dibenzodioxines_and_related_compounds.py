"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds

Checks if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl, or related compound by identifying core structures and counting chlorine or bromine substituents.

This improved version addresses previous issues by using more general SMARTS patterns for core structures
and implementing size constraints to exclude large molecules not relevant to the class.
"""

from rdkit import Chem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl,
    or related compound based on its SMILES string.

    These compounds are characterized by the presence of
    dibenzo-p-dioxin, dibenzofuran, or biphenyl core structures, and multiple chlorine or bromine atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is in the class, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Limit the size of the molecule to exclude large, complex molecules
    max_heavy_atoms = 50
    num_heavy_atoms = mol.GetNumHeavyAtoms()
    if num_heavy_atoms > max_heavy_atoms:
        return False, f"Molecule has {num_heavy_atoms} heavy atoms; exceeds limit of {max_heavy_atoms}"

    # Limit the number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 4:
        return False, f"Molecule has {num_rings} rings; exceeds limit of 4"

    # Define core structure patterns using SMARTS
    # Use more general patterns for better matching
    dibenzo_dioxin_pattern = Chem.MolFromSmarts('c1cc2oc3ccccc3oc2cc1')  # Dibenzo-p-dioxin core
    dibenzofuran_pattern = Chem.MolFromSmarts('c1cc2oc3ccccc3c2cc1')    # Dibenzofuran core
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c2ccccc2')          # Biphenyl core

    # Check for core structures
    has_dibenzo_dioxin = mol.HasSubstructMatch(dibenzo_dioxin_pattern)
    has_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)

    if has_dibenzo_dioxin:
        core_structure = 'dibenzo-p-dioxin'
    elif has_dibenzofuran:
        core_structure = 'dibenzofuran'
    elif has_biphenyl:
        core_structure = 'biphenyl'
    else:
        return False, "Molecule does not contain dibenzo-p-dioxin, dibenzofuran, or biphenyl core structures"

    # Count number of chlorine and bromine atoms
    num_cl = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 17)
    num_br = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 35)
    total_halogens = num_cl + num_br

    if total_halogens < 2:
        return False, f"Molecule contains {total_halogens} halogen atom(s); requires at least 2 for classification"

    # Ensure that chlorine and bromine are the only halogens present
    num_f = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
    num_i = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 53)
    if num_f > 0 or num_i > 0:
        return False, "Molecule contains fluorine or iodine atoms; only chlorine and bromine are allowed"

    # Ensure that the molecule does not contain other large functional groups by checking for other heteroatoms
    # For simplicity, we can limit heteroatoms to oxygen, chlorine, and bromine
    allowed_atomic_numbers = {6, 1, 8, 17, 35}  # Carbon, hydrogen, oxygen, chlorine, bromine
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_numbers:
            return False, f"Molecule contains atom {atom.GetSymbol()}; only C, H, O, Cl, and Br are allowed"

    return True, f"Molecule contains {core_structure} core and is polychlorinated/polybrominated with {total_halogens} halogen atoms"