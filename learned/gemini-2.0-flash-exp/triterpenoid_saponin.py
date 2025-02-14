"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
"""
Classifies: triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin based on its SMILES string.
    A triterpenoid saponin is a terpene glycoside in which the terpene moiety is a triterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for triterpenoid core (fused or bridged 6-membered rings)
    # A flexible pattern to match fused or bridged 6 member rings
    triterpenoid_core_pattern = Chem.MolFromSmarts("[CX4]1[CX4][CX4][CX4][CX4][CX4]1[CX4]~[CX4]2[CX4][CX4][CX4][CX4][CX4]2")
    core_matches = mol.GetSubstructMatches(triterpenoid_core_pattern)
    if not core_matches:
       return False, "No triterpenoid core detected"


    # Count carbons and oxygens in triterpenoid core
    core_match_atoms = [atom for match in core_matches for atom_idx in match for atom in mol.GetAtoms() if atom.GetIdx() == atom_idx]
    core_c_count = sum(1 for atom in core_match_atoms if atom.GetAtomicNum() == 6)
    core_o_count = sum(1 for atom in core_match_atoms if atom.GetAtomicNum() == 8)
    if core_c_count < 25: # minimum carbon atoms for a triterpenoid core
        return False, "Triterpenoid core has too few carbons"


    # Look for glycosidic bonds connecting triterpenoid core or ring system to a carbohydrate (pyranose/furanose)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[CX4;R]~[OX2]~[CX4;R1]")  #C-O-C where both C are in a ring
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    if len(glycosidic_matches) < 1:
        return False, "No glycosidic bond detected"


    # Check if attached groups to glycosidic bonds are sugars (simplified check using ring size and O count)
    sugar_pattern1 = Chem.MolFromSmarts("[OX2]1[CX4][CX4][CX4][CX4][CX4]1") # 6-membered ring sugar
    sugar_pattern2 = Chem.MolFromSmarts("[OX2]1[CX4][CX4][CX4][CX4]1")  #5-membered ring sugar
    sugar_found = False
    for match in glycosidic_matches:
      for atom_index in match:
          atom = mol.GetAtomWithIdx(atom_index)
          if atom.GetAtomicNum() == 8:
              for neighbor in atom.GetNeighbors():
                  if neighbor.GetIdx() == match[0] or neighbor.GetIdx() == match[2]: # don't consider atoms of the glycosidic bond
                     continue
                  submol = Chem.PathToSubmol(mol, [atom_index, neighbor.GetIdx()])
                  if submol.HasSubstructMatch(sugar_pattern1) or submol.HasSubstructMatch(sugar_pattern2) :
                     sugar_found = True
                     break
      if sugar_found:
          break

    if not sugar_found:
      return False, "Glycosidic bond does not appear to be connected to a sugar moiety"

    return True, "Contains a triterpenoid core with glycosidic bond(s) and sugar(s)."