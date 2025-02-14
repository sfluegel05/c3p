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

    # Core structure patterns, allowing for substitutions on the rings
    # Dibenzodioxin core
    dioxin_core_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3]c1-O-[cX3]2[cX3][cX3][cX3][cX3]c2-O-")
    # Dibenzofuran core
    furan_core_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3]c1-O-[cX3]2[cX3][cX3][cX3][cX3]c2")
    # Biphenyl core
    biphenyl_core_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3]c1-[cX3]2[cX3][cX3][cX3][cX3]c2")


    # Check for core structures
    if mol.HasSubstructMatch(dioxin_core_pattern):
      # Count aromatic halogens (Cl or Br)
      halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35] and atom.GetIsAromatic())
      if halogen_count >= 2:
          return True, "Polychlorinated dibenzodioxin detected"
      else:
          return False, "Dibenzodioxin core found, but lacks sufficient aromatic halogens"

    if mol.HasSubstructMatch(furan_core_pattern):
      # Count aromatic halogens (Cl or Br)
      halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35] and atom.GetIsAromatic())
      if halogen_count >= 2:
          return True, "Polychlorinated dibenzofuran detected"
      else:
        return False, "Dibenzofuran core found, but lacks sufficient aromatic halogens"

    if mol.HasSubstructMatch(biphenyl_core_pattern):
      # Count aromatic halogens (Cl or Br)
      halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35] and atom.GetIsAromatic())
      if halogen_count >= 2:
          return True, "Polychlorinated or polybrominated biphenyl detected"
      else:
          return False, "Biphenyl core found, but lacks sufficient aromatic halogens"

    #Check for biphenyl core with aromatic halogens (in case the core is not a simple 2 rings like in the previous substructure).
    any_biphenyl_pattern = Chem.MolFromSmarts("[cX3]~[cX3]")
    if mol.HasSubstructMatch(any_biphenyl_pattern):
        # Count aromatic halogens (Cl or Br)
        halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35] and atom.GetIsAromatic())
        if halogen_count >= 2:
          return True, "Biphenyl with aromatic halogens detected"

    return False, "Does not match any of the defined substructures or it lacks sufficient aromatic halogens"