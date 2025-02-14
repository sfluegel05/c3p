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

    # 1. Define a flexible steroid core pattern with a marked position (position 3)
    # This pattern matches the four fused rings of a steroid, with possible double bonds
    # and marks a carbon at position 3 using the special atom index 1
    steroid_core_pattern = Chem.MolFromSmarts("[C]1-[C]-[C]-[C]2-[C]-[C]-[C]3-[C]-[C]4-[C]1-[C]243")

    # Find the steroid core in the molecule
    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if not core_match:
        return False, "No steroid core found"
    
    #get position 3, by the marked atom in the SMARTS
    pos3_atom = core_match[0]

    # 2. Define the ester pattern connected to any C of the steroid core
    ester_pattern = Chem.MolFromSmarts("[C]1-[C]-[C]-[C]2-[C]-[C]-[C]3-[C]-[C]4-[C]1-[C]243[C](=O)O")

    # Find all ester groups
    ester_matches = mol.GetSubstructMatch(ester_pattern)


    if not ester_matches:
      return False, "No ester group found connected to the steroid core"
    
    return True, "Sterol ester detected: Contains a steroid core with an ester bond at position 3"