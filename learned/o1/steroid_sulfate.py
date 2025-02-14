"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
"""

from rdkit import Chem
from rdkit.Chem import rdqueries

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a steroid molecule where a hydroxy group of the steroid core
    is esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid core pattern (cyclopenta[a]phenanthrene skeleton)
    steroid_smarts = '[#6]12[#6][#6][#6]3[#6]([#6][#6][#6]([#6]4[#6][#6][#6]([#6]1)[#6]([#6]4)[#6]([#6]2)[#6]3)[#6]([#6]([#6][#6][#6][#6])[#6][#6])'
    steroid_core = Chem.MolFromSmarts(steroid_smarts)

    # Check for steroid core
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core found"

    # Define sulfate ester group pattern (attached via oxygen)
    sulfate_ester_smarts = '[O;D2]-S(=O)(=O)-[O;D1]'
    sulfate_ester = Chem.MolFromSmarts(sulfate_ester_smarts)

    # Find sulfate ester groups
    sulfate_matches = mol.GetSubstructMatches(sulfate_ester)
    if not sulfate_matches:
        return False, "No sulfate ester group found"

    # Check that sulfate group is connected to steroid core via oxygen
    steroid_match_atoms = mol.GetSubstructMatch(steroid_core)
    steroid_atom_indices = set(steroid_match_atoms)

    for match in sulfate_matches:
        # The first oxygen in the pattern is connected to the steroid
        ester_oxygen_idx = match[0]
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)

        connected = False
        for neighbor in ester_oxygen_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in steroid_atom_indices:
                connected = True
                break

        if connected:
            return True, "Contains steroid core with sulfate ester group attached via oxygen"

    return False, "Sulfate group not attached to steroid core via ester linkage"