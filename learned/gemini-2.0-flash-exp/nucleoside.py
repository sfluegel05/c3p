"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside consists of a nucleobase and a ribose or deoxyribose sugar linked by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ribose/deoxyribose sugar pattern (5-membered ring with 1 oxygen)
    sugar_pattern = Chem.MolFromSmarts("[CX4]1[OX2][CX4][CX4][CX4]1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No ribose or deoxyribose sugar found."

    # Nucleobase patterns (purine or pyrimidine) - more general patterns
    purine_pattern = Chem.MolFromSmarts("[n]1cncn2c1cc[n]2")
    pyrimidine_pattern = Chem.MolFromSmarts("[n]1cc[n]c(=[OX1])c1")

    purine_matches = mol.GetSubstructMatches(purine_pattern)
    pyrimidine_matches = mol.GetSubstructMatches(pyrimidine_pattern)

    if not purine_matches and not pyrimidine_matches:
      return False, "No purine or pyrimidine nucleobase found."

    # Check for glycosidic bond (C-N bond linking sugar and base).
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4]-[NX3]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if not glycosidic_matches:
       return False, "No C-N glycosidic bond found."
    
    #Check that the glycosidic bond carbon is part of the sugar
    sugar_atoms = set()
    for match in sugar_matches:
        for a_idx in match:
            sugar_atoms.add(a_idx)

    connected = False
    for match in glycosidic_matches:
      glycosidic_c = match[0]
      if glycosidic_c in sugar_atoms:
          connected = True
          break
    if not connected:
        return False, "Sugar and base not connected via glycosidic bond"


    return True, "Contains a nucleobase linked to a ribose or deoxyribose sugar."