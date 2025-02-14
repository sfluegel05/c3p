"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
"""

from rdkit import Chem

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
    steroid_core_smarts = '[#6]1~[#6]~[#6]2~[#6](~[#6]~[#6]~[#6]~2)~[#6]3~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6]1'
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)

    # Check for steroid core
    steroid_matches = mol.GetSubstructMatches(steroid_core)
    if not steroid_matches:
        return False, "No steroid core found"

    # Define sulfate ester group pattern (attached via oxygen)
    sulfate_ester_smarts = 'OS(=O)(=O)[O]'
    sulfate_ester = Chem.MolFromSmarts(sulfate_ester_smarts)

    # Find sulfate ester groups
    sulfate_matches = mol.GetSubstructMatches(sulfate_ester)
    if not sulfate_matches:
        return False, "No sulfate ester group found"

    # For each sulfate group, check if it's connected to the steroid core via oxygen
    for sulfate_match in sulfate_matches:
        ester_oxygen_idx = sulfate_match[0]  # Oxygen attached to sulfur
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
        connected_to_steroid = False

        for neighbor in ester_oxygen_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            # Check if neighbor atom is part of steroid core
            for steroid_match in steroid_matches:
                if neighbor_idx in steroid_match:
                    connected_to_steroid = True
                    break
            if connected_to_steroid:
                break

        if connected_to_steroid:
            return True, "Contains steroid core with sulfate ester group attached via oxygen"

    return False, "Sulfate group not attached to steroid core via ester linkage"