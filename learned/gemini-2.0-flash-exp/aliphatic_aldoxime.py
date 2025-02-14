"""
Classifies: CHEBI:82744 aliphatic aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aliphatic_aldoxime(smiles: str):
    """
    Determines if a molecule is an aliphatic aldoxime based on its SMILES string.
    An aliphatic aldoxime is an aldoxime where the carbon attached to the oxime group is derived from an aliphatic aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic aldoxime, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check for the sp2 oxime group (C=N-O) and the directly linked carbon must be aldehyde-like
    oxime_pattern = Chem.MolFromSmarts("[CX3H1]=[NX2]-O")
    matches = mol.GetSubstructMatches(oxime_pattern)
    if not matches:
      return False, "No sp2 aldoxime group found"

    for match in matches:
      carbon_index = match[0]
      carbon_atom = mol.GetAtomWithIdx(carbon_index)
      is_aliphatic = True
    
      for neighbor in carbon_atom.GetNeighbors():
          if neighbor.GetAtomicNum() != 1 and neighbor.GetAtomicNum() != 6 and neighbor.GetIdx() != match[1]:  #check that it is either C or H, or the nitrogen of the oxime.
              is_aliphatic = False
              break
      if not is_aliphatic:
          return False, "Aldoxime carbon is not directly connected to an aliphatic carbon"


    # 4. Return True if all conditions are met
    return True, "Molecule is an aliphatic aldoxime"