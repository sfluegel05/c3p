"""
Classifies: CHEBI:134045 polychlorinated dibenzodioxines and related compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polychlorinated_dibenzodioxines_and_related_compounds(smiles: str):
    """
    Determines if a molecule is a polychlorinated dibenzodioxine, dibenzofuran, or polychlorinated/polybrominated biphenyl, or a related analog.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule belongs to the class, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count all halogens (Cl or Br)
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35])
    if halogen_count < 2:
        return False, "Less than 2 halogens found"

    # Core structure patterns, allowing for substitutions on the rings
    # Dibenzodioxin core
    dioxin_core_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c]c1-O-[c]2[c][c][c][c]c2-O-")
    # Dibenzofuran core
    furan_core_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c]c1-O-[c]2[c][c][c][c]c2")
    # Biphenyl core (any two aromatic rings connected)
    biphenyl_core_pattern = Chem.MolFromSmarts("[c]~[c]")
    # Ether biaryl
    ether_pattern = Chem.MolFromSmarts("[c]~[O]~[c]")

    # Check for core structures
    if mol.HasSubstructMatch(dioxin_core_pattern):
        return True, "Polychlorinated dibenzodioxin detected"

    if mol.HasSubstructMatch(furan_core_pattern):
        return True, "Polychlorinated dibenzofuran detected"

    if mol.HasSubstructMatch(biphenyl_core_pattern):
      return True, "Polychlorinated or polybrominated biphenyl detected"

    if mol.HasSubstructMatch(ether_pattern):
      return True, "Polychlorinated biaryl ether detected"

    if mol.HasSubstructMatch(biphenyl_core_pattern) and mol.HasSubstructMatch(ether_pattern):
        return True, "Polychlorinated biphenyl ether"

    return False, "Does not match any of the defined substructures or it lacks sufficient halogens"