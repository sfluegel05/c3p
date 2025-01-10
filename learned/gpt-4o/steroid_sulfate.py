"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester of a steroid, featuring a steroid backbone
    with an attached sulfate group.

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

    # Pattern for sulfuric ester group (e.g., various sulfate representations)
    sulfate_group_patterns = [
        Chem.MolFromSmarts("OS(=O)(=O)O"),
        Chem.MolFromSmarts("OS(=O)(=O)[O-]"),
        Chem.MolFromSmarts("[O-]S(=O)(=O)O")
    ]
    
    # Check for the presence of any recognizable sulfate group pattern
    sulfate_found = any(mol.HasSubstructMatch(pattern) for pattern in sulfate_group_patterns)
    if not sulfate_found:
        return False, "No recognizable sulfate group found"
    
    # Pattern for detecting a more flexible steroid core
    steroid_core_pattern = Chem.MolFromSmarts("C1C2CC3CC(C1)CCC3C4=C2C=CC=C4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core structure found"

    # Confirm at least one sulfate group is attached correctly
    attachment_found = False
    for pattern in sulfate_group_patterns:
        sulfur_matches = mol.GetSubstructMatches(pattern)
        for match in sulfur_matches:
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                # Check for attachment to steroid-like carbon atom (part of the core)
                for neighbor in atom.GetNeighbors():
                    if neighbor.IsInRing() and neighbor.GetAtomicNum() == 6:
                        # Confirm it is within the steroid core
                        if any(neighbor.HasSubstructMatch(steroid_core_pattern) for neighbor in atom.GetNeighbors()):
                            attachment_found = True
                            break
                if attachment_found:
                    break
            if attachment_found:
                break
    
    if not attachment_found:
        return False, "Sulfate group not correctly attached to steroid structure"
    
    return True, "Contains steroid structure with sulfate group attached correctly"

# Example usage (test with one of the SMILES from the list above):
print(is_steroid_sulfate("C1=C(OS(O)(=O)=O)C=CC2=C1CC[C@]3([C@@]4(CC[C@]([C@]4(CC[C@@]32[H])C)(C#C)O)[H])[H]"))