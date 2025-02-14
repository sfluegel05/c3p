"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid has the basic steroid core and a beta-configured hydroxyl group at position 17.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the steroid core (more flexible)
    # using atom mapping to identify position C17, labeled as 17
    steroid_core_pattern = Chem.MolFromSmarts("[C:1]12[C:2]([C:3][C:4]([C:1])([C:5]([C:6]([C:7]([C:2])[C:8]([C:3])C)CC[C:4])C)C")
    match = mol.GetSubstructMatch(steroid_core_pattern)
    if not match:
        return False, "Molecule does not have the basic steroid core structure."

    # Get the C17 atom
    c17_idx = match[7]
    c17_atom = mol.GetAtomWithIdx(c17_idx)
    
    # Check for exactly one hydroxyl neighbor on C17
    oh_count = 0
    oxygen_neighbor = None
    for nbr in c17_atom.GetNeighbors():
      if nbr.GetAtomicNum() == 8:
        oh_count += 1
        oxygen_neighbor = nbr

    if oh_count != 1:
        return False, "C17 does not have exactly one hydroxyl group"
    
    # Check the stereochemistry of the C17 atom
    if c17_atom.HasProp('_CIPCode'):
        stereo = c17_atom.GetProp('_CIPCode')
    elif c17_atom.HasProp('_Chirality'):
        stereo = c17_atom.GetProp('_Chirality')
    else:
        return False, "C17 has no stereochemistry defined."
        
    # Check if the stereochemistry of carbon is R
    if stereo != "R":
      return False, "The hydroxyl group at position 17 is not beta-configured"
    
    
    return True, "Molecule contains a steroid core with a beta-hydroxy group at position 17."