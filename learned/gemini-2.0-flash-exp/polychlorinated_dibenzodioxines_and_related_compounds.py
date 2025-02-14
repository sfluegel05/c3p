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

    # Count aromatic halogens (Cl or Br)
    halogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [17, 35] and atom.GetIsAromatic())
    if halogen_count < 2:
        return False, "Must have at least two aromatic chlorine or bromine atoms."


    # Check for dibenzodioxin, dibenzofuran or biphenyl core structure, allowing for substitutions.
    dioxin_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3][cX3]1-O-[cX3]1[cX3][cX3][cX3][cX3][cX3]1-O-[cX3]1[cX3][cX3][cX3][cX3][cX3]1")
    furan_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3][cX3]1-O-[cX3]1[cX3][cX3][cX3][cX3][cX3]1-[cX3]1[cX3][cX3][cX3][cX3][cX3]1")
    biphenyl_pattern = Chem.MolFromSmarts("[cX3]1[cX3][cX3][cX3][cX3][cX3]1-[cX3]1[cX3][cX3][cX3][cX3][cX3]1")

    if mol.HasSubstructMatch(dioxin_pattern):
      return True, "Polychlorinated dibenzodioxin detected"
    if mol.HasSubstructMatch(furan_pattern):
        return True, "Polychlorinated dibenzofuran detected"
    if mol.HasSubstructMatch(biphenyl_pattern):
      return True, "Polychlorinated or polybrominated biphenyl detected"

    #Check for any biphenyl core with aromatic halogens
    any_biphenyl_pattern = Chem.MolFromSmarts("[cX3]~[cX3]")
    if mol.HasSubstructMatch(any_biphenyl_pattern) and halogen_count >= 2:
      return True, "Biphenyl with aromatic halogens detected"
      
    return False, "Does not match any of the defined substructures or it lacks sufficient aromatic halogens"