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

    # Expand potential patterns for a steroid backbone
    steroid_patterns = [
        Chem.MolFromSmarts("[C@@H]1CC[C@@]2[C@H]3CC[C@]4(C)[C@]3(CC[C@H]2[C@@H]1C)O4"),  # More detailed steroid core
        Chem.MolFromSmarts("[C@H]1CC[C@@]2(C[C@@H](CC[C@@]2C1)O)C"),  # Core with oxygen/hydroxyl groups included
    ]

    # Check if there's any pattern matching a steroid backbone
    steroid_match = any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns)
    if not steroid_match:
        return False, "No steroid backbone found"

    # SMARTS pattern to identify sulfate group bound through an ester linkage (with slight flexibility)
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")

    # Look for sulfate groups connected to hydroxyl groups
    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern)
    for match in sulfate_matches:
        oxygen_atom = [mol.GetAtomWithIdx(idx) for idx in match if mol.GetAtomWithIdx(idx).GetSymbol() == 'O'][0]
        carbon_neighbors = [nbr for nbr in oxygen_atom.GetNeighbors() if nbr.GetSymbol() == 'C']
        
        for carbon in carbon_neighbors:
            # Check if carbon originally belongs to hydroxy group of steroid (OH pattern)
            if carbon.GetDegree() == 3:  # Simple way to ensure it's part of a larger chain, not dangling
                for o_neighbor in carbon.GetNeighbors():
                    if o_neighbor.GetSymbol() == 'O' and o_neighbor != oxygen_atom:
                        if any(nbr.GetSymbol() == 'H' for nbr in o_neighbor.GetNeighbors()):
                            return True, "Contains sulfate ester linked to steroid backbone at a hydroxy group"

    return False, "Sulfate groups found but not linked properly to steroid backbone"

# Example for testing
smiles = "[Na+].[H][C@]12CC[C@]3(C)C(=O)CC[C@@]3([H])[C@]1([H])CCc1cc(OS([O-])(=O)=O)ccc21"
result, reason = is_steroid_sulfate(smiles)
print(f"Result: {result}, Reason: {reason}")