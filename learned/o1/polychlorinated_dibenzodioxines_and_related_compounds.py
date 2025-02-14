"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: polychlorinated dibenzodioxines and related compounds

Checks if a molecule is a polychlorinated dibenzodioxin, dibenzofuran, biphenyl, or related compound by identifying core structures and counting halogen substituents.
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

    # Define core structure patterns using SMARTS
    dibenzo_dioxin_pattern = Chem.MolFromSmarts('O1c2ccccc2Oc3ccccc13')  # Dibenzo-p-dioxin core
    dibenzofuran_pattern = Chem.MolFromSmarts('O1c2ccccc2c3ccccc13')    # Dibenzofuran core
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

    return True, f"Molecule contains {core_structure} core and is polychlorinated/polybrominated with {total_halogens} halogen atoms"