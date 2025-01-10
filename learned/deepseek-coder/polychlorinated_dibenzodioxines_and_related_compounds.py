"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
"""
Classifies: Polychlorinated dibenzodioxins and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxin or related compound based on its SMILES string.
    These compounds include polychlorinated dibenzodioxins, dibenzofurans, and biphenyls with multiple halogen substitutions.

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

    # Define core structures for dibenzodioxin, dibenzofuran, and biphenyl
    dibenzodioxin_pattern = Chem.MolFromSmarts("O1c2ccccc2c3ccccc3O1")
    dibenzofuran_pattern = Chem.MolFromSmarts("O1c2ccccc2c3ccccc31")
    biphenyl_pattern = Chem.MolFromSmarts("c1ccccc1-c2ccccc2")

    # Check if the molecule contains any of the core structures
    has_dibenzodioxin = mol.HasSubstructMatch(dibenzodioxin_pattern)
    has_dibenzofuran = mol.HasSubstructMatch(dibenzofuran_pattern)
    has_biphenyl = mol.HasSubstructMatch(biphenyl_pattern)

    if not (has_dibenzodioxin or has_dibenzofuran or has_biphenyl):
        return False, "No dibenzodioxin, dibenzofuran, or biphenyl core structure found"

    # Count the number of halogen atoms (Cl or Br)
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35])
    if halogen_count < 2:
        return False, f"Found {halogen_count} halogen atoms, need at least 2"

    # Check for multiple halogen substitutions (at least 2 halogens on the core structure)
    if has_dibenzodioxin:
        core_matches = mol.GetSubstructMatches(dibenzodioxin_pattern)
    elif has_dibenzofuran:
        core_matches = mol.GetSubstructMatches(dibenzofuran_pattern)
    else:
        core_matches = mol.GetSubstructMatches(biphenyl_pattern)

    core_halogen_count = 0
    for match in core_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() in [17, 35]:
                core_halogen_count += 1

    if core_halogen_count < 2:
        return False, f"Found {core_halogen_count} halogen atoms on the core structure, need at least 2"

    return True, "Contains a dibenzodioxin, dibenzofuran, or biphenyl core with multiple halogen substitutions"