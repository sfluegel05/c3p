"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is a steroid with a hydroxyl group at the 3-beta position.
    It has a tetracyclic ring system closely related to cholestan-3-ol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid core SMARTS pattern (tetracyclic ring system with 6-6-6-5 fusion), includes numbering of key carbons
    steroid_core_pattern = Chem.MolFromSmarts("[C]12[C]([C]([C]([C]([C]1)[C]3)C)CC[C]4)[C]([C](C2)CC[C]5[C]3([C](C([C]5)C)C))[C]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"

    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if not core_match:
        return False, "No match of steroid_core_pattern"


    # Identify C3 atom (index 2 in the match since it is the second carbon in the SMARTS)
    c3 = mol.GetAtomWithIdx(core_match[2])

    # Check for 3-beta hydroxyl group. Must be attached to c3 and have one H.
    has_3beta_oh = False
    for neighbor in c3.GetNeighbors():
        if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
            has_3beta_oh = True
            break
    
    if not has_3beta_oh:
        return False, "No hydroxyl group at the 3-beta position"

    # Identify C17 atom (index 17 in the match)
    c17 = mol.GetAtomWithIdx(core_match[17])

    # Check for side chain at C17 (at least one carbon atom attached)
    has_sidechain = False
    for neighbor in c17.GetNeighbors():
        if neighbor.GetSymbol() == 'C':
           has_sidechain = True
           break
    
    if not has_sidechain:
      return False, "No side chain at C17"


    return True, "Contains a steroid core, a hydroxyl group at the 3-beta position, and a side chain at C17."