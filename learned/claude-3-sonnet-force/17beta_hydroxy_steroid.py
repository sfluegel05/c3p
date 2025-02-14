"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: CHEBI:35577 17beta-hydroxy steroid
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is a steroid with a hydroxy group at position 17 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone
    steroid_pattern = Chem.MolFromSmarts("[C&r1,r2,r3,r4]") # Any atom in a ring system with 4 rings
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Find the 17-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@](O)[H]") # 17-hydroxy group in beta configuration
    matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not matches:
        return False, "No 17-hydroxy group found"

    # Check if the 17-hydroxy group is part of the steroid backbone
    ring_info = mol.GetRingInfo()
    for match in matches:
        atom = mol.GetAtomWithIdx(match[0])
        if any(ring_info.IsBondInRingOfSize(bond.GetIdx(), 5) or ring_info.IsBondInRingOfSize(bond.GetIdx(), 6) for bond in atom.GetBonds()):
            return True, "Molecule has a 17beta-hydroxy group on the steroid backbone"

    return False, "17-hydroxy group is not part of the steroid backbone"