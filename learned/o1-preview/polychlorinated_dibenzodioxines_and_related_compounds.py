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
    and have at least two chlorine or bromine atoms.

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

    # Define core structures using SMILES
    biphenyl_smiles = 'c1ccccc1-c2ccccc2'  # Biphenyl core
    dibenzodioxin_smiles = 'O1c2ccccc2Oc3ccccc13'  # Dibenzodioxin core
    dibenzofuran_smiles = 'c1ccc2c(c1)oc1ccccc12'  # Dibenzofuran core

    biphenyl_mol = Chem.MolFromSmiles(biphenyl_smiles)
    dibenzodioxin_mol = Chem.MolFromSmiles(dibenzodioxin_smiles)
    dibenzofuran_mol = Chem.MolFromSmiles(dibenzofuran_smiles)

    # Check for core structures
    core_matches = []
    if mol.HasSubstructMatch(biphenyl_mol):
        core_matches.append('biphenyl')
    if mol.HasSubstructMatch(dibenzodioxin_mol):
        core_matches.append('dibenzodioxin')
    if mol.HasSubstructMatch(dibenzofuran_mol):
        core_matches.append('dibenzofuran')

    if not core_matches:
        return False, "Does not contain dibenzodioxin, dibenzofuran, or biphenyl core structure"

    # Count number of Cl and Br atoms in the molecule
    num_cl = 0
    num_br = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Cl':
            num_cl +=1
        elif atom.GetSymbol() == 'Br':
            num_br +=1

    num_halogen = num_cl + num_br

    if num_halogen < 2:
        return False, f"Contains core structure(s) {', '.join(core_matches)}, but has less than 2 Cl or Br atoms"

    return True, f"Contains {', '.join(core_matches)} core with {num_halogen} Cl/Br atoms"