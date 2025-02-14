"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid with an ester bond at the 3-hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Define a steroid core pattern with a marked position (position 3)
    # This pattern matches the four fused rings of a steroid, with possible double bonds
    # and it marks the carbon at position 3 using the special atom index 1
    steroid_core_pattern = Chem.MolFromSmarts("[CX3,CX4]1[CX3,CX4][CX3,CX4][CX3,CX4]2[CX3,CX4][CX3,CX4][CX3,CX4]3[CX3,CX4][CX3,CX4]4[CX3,CX4]1[CX3,CX4]243")


    # Find the steroid core in the molecule
    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if not core_match:
        return False, "No steroid core found"
    
    #get position 3, by the marked atom in the SMARTS
    pos3_atom = core_match[0]

    # 2. Define the ester pattern
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    
    # Find all ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # 3. Check if any ester is connected to position 3
    ester_at_pos3 = False
    for ester_match in ester_matches:
        for ester_atom in ester_match:
          ester_atom_obj = mol.GetAtomWithIdx(ester_atom)
          for neighbor_atom in ester_atom_obj.GetNeighbors():
              if neighbor_atom.GetIdx() == pos3_atom:
                ester_at_pos3 = True
                break
        if ester_at_pos3:
          break
    
    if not ester_at_pos3:
      return False, "No ester group found at position 3 of the steroid core"


    return True, "Sterol ester detected: Contains a steroid core with an ester bond at position 3"