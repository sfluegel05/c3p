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

    # Steroid core SMARTS pattern (tetracyclic ring system)
    steroid_core_pattern = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]~[C]~[C]12~[C]~[C]~[C]~[C]~[C]23~[C]~[C]~[C]~[C]34~[C]~[C]~[C]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found"
    
    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    if not core_match:
        return False, "No match of steroid_core_pattern"

    # Find the 3-beta carbon (the one with the OH group)
    has_3beta_oh = False
    c17 = None
    for atom_index in core_match:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.IsInRingSize(6) and atom.IsInRingSize(6) and atom.IsInRingSize(5):
                for neighbor in atom.GetNeighbors():
                     if neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1:
                        has_3beta_oh = True
                        break
                if has_3beta_oh:
                    break

    if not has_3beta_oh:
        return False, "No hydroxyl group at the 3-beta position"

    # Find C17 (carbon connected to the 5-membered ring, will have a side chain)
    c17_candidates = []
    for atom_index in core_match:
         atom = mol.GetAtomWithIdx(atom_index)
         if atom.IsInRingSize(5) :
             c17_candidates.append(atom)
    
    has_sidechain = False

    if len(c17_candidates) == 1: # This is the correct approach if one 5 membered ring carbon is involved in ring fusion.
         c17 = c17_candidates[0]
    else:
      return False, "Unable to identify C17 carbon"

    for neighbor in c17.GetNeighbors():
        if neighbor.GetSymbol() == 'C' and neighbor.GetTotalNumHs() <= 2: # can be CH, CH2 but not CH3
            has_sidechain = True
            break
            
    if not has_sidechain:
      return False, "No side chain at C17"

    return True, "Contains a steroid core, a hydroxyl group at the 3-beta position, and a side chain at C17."