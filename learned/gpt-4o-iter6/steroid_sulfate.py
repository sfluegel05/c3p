"""
Classifies: CHEBI:16158 steroid sulfate
"""
from rdkit import Chem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a sulfuric ester obtained by the condensation of a hydroxy group of any steroid with sulfuric acid.

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

    # Basic steroid core patterns, deducing a typical fused-ring backbone system
    core_patterns = [
        "C1CCC2C(C1)CCC3C2CCC4C3CCC4",  # Typical backbone
        "[C@H]1([C@H]2CC[C@@H]3[C@@H]2[C@@H](CC[C@@H]4[C@@H]3CC[C@@H]1[C@@H]4CO)O)"
        # These can be improved/extended based on specific examples or libraries
    ]

    # Check presence of a steroid-like core
    steroid_match = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in core_patterns)
    if not steroid_match:
        return False, "No steroid backbone found"

    # Recognize sulfate group
    sulfate_patterns = [
        "OS(=O)(=O)[O-]",  # Handling of sulfate ions
        "OS([O-])(=O)=O",  # Additional potential form in salts
        "OS(=O)(=O)O"     # Sulfate ester
    ]

    # Check any pattern of sulfate ester presence
    sulfate_match = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in sulfate_patterns)
    if not sulfate_match:
        return False, "No sulfate groups correctly attached"

    # Ensuring sulfate bond attach to potential hydroxy site's carbon first before linkage
    alcohol_attached_sulfate = False
    for sulfate_pattern in sulfate_patterns:
        sulfate = Chem.MolFromSmarts(sulfate_pattern)
        matches = mol.GetSubstructMatches(sulfate)
        for match in matches:
            # Should possess connectivity to alcohol site that forms the ester
            atom = mol.GetAtomWithIdx(match[0])  # Getting S atom
            neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O' and n.GetDegree() == 2]  # Check O attaches to C
            for n in neighbors:
                if any(nei.GetSymbol() == 'C' for nei in n.GetNeighbors()):
                    alcohol_attached_sulfate = True
                    break

    # If all checks succeed, return true.
    if alcohol_attached_sulfate:
        return True, "Contains sulfate ester linked to steroid backbone at a hydroxy group"

    return False, "Sulfate groups found but not linked properly to steroid backbone"

# Example for testing
smiles = "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O"
result, reason = is_steroid_sulfate(smiles)
print(f"Result: {result}, Reason: {reason}")