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

    # Expanded steroid core patterns - including variations found commonly in steroids
    core_patterns = [
        "C1CCC2C(C1)CCC3C2CCC4C3CCC4",  # Basic steroid backbone
        "[C@H]1([C@@H]2CCC3C=C(C2)CC4[C@@H]1CC(O)CC34)",  # Additional possible steroid core
        "C1C2CC3CCC4C(CCC4C3)C2C1"  # Simplified pattern with stereochemistry consideration
    ]

    # Check presence of a steroid-like core
    steroid_match = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in core_patterns)
    if not steroid_match:
        return False, "No steroid backbone found"

    # Sulfate group matching, focusing on ester linkages
    sulfate_patterns = [
        "O[S](=O)(=O)[O-]",   # With negative charge for salts
        "O[S](=O)(=O)O",      # Neutral sulfate ester
        "O[S](=O)(=O)[O0]"    # Consider neutral sulfate without explicit charge
    ]

    # Verification of sulfate connected to hydroxy group
    sulfate_found = False
    for sulfate_pattern in sulfate_patterns:
        sulfate = Chem.MolFromSmarts(sulfate_pattern)
        matches = mol.GetSubstructMatches(sulfate)
        sulfate_group_atoms = [mol.GetAtomWithIdx(match[0]) for match in matches]
        
        for sulfur_atom in sulfate_group_atoms:
            oxygen_neighbors = [n for n in sulfur_atom.GetNeighbors() if n.GetSymbol() == 'O']
            for oxygen in oxygen_neighbors:
                carbon_neighbors = [n for n in oxygen.GetNeighbors() if n.GetSymbol() == 'C']
                for carbon in carbon_neighbors:
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("CO")):  # Check if connected to OH group
                        for neighbor in carbon.GetNeighbors():
                            if neighbor.GetSymbol() == 'O' and 'H' in [a.GetSymbol() for a in neighbor.GetNeighbors()]:
                                sulfate_found = True
                                break

    if not sulfate_found:
        return False, "Sulfate groups found but not linked properly to steroid backbone"

    return True, "Contains sulfate ester linked to steroid backbone at a hydroxy group"

# Example for testing
smiles = "[Na+].C[C@]12CC[C@H]3C(=CCc4cc(OS([O-])(=O)=O)ccc34)[C@@H]1CCC2=O"
result, reason = is_steroid_sulfate(smiles)
print(f"Result: {result}, Reason: {reason}")