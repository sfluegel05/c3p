"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for monosaccharide (cyclic or acyclic)
    # This pattern matches sugars with multiple hydroxyl groups
    sugar_patterns = [
        "[C;!R][C;!R][C;!R][C;!R][C;!R][C;!R]",  # Linear hexose sugar
        "[C;R][C;R][C;R][C;R][O;R]",             # 5-membered cyclic sugar
        "[C;R][C;R][C;R][C;R][C;R][O;R]",        # 6-membered cyclic sugar
    ]
    sugar_found = False
    for smarts in sugar_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            sugar_found = True
            break

    if not sugar_found:
        return False, "No monosaccharide unit found"

    # Define SMARTS pattern for phosphate ester group attached via oxygen
    phosphate_patterns = [
        "[O]-P(=O)([O])[O]",    # Phosphate ester linkage
        "[O]-P(=O)([O-])[O-]",  # Phosphate ester with negative charges
        "[O]-P(=O)([O])[O-]",
    ]
    phosphate_found = False
    for smarts in phosphate_patterns:
        phosphate_pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(phosphate_pattern):
            phosphate_found = True
            break

    if not phosphate_found:
        return False, "No phosphate ester group found"

    # Check that the phosphate group is attached to the sugar via an oxygen atom
    # Find the indices of the sugar atoms
    sugar_atoms = set()
    for smarts in sugar_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            sugar_atoms.update(match)

    # Find the phosphate ester groups
    phosphate_atoms = set()
    for smarts in phosphate_patterns:
        phosphate_pattern = Chem.MolFromSmarts(smarts)
        matches = mol.GetSubstructMatches(phosphate_pattern)
        for match in matches:
            phosphate_atoms.add(match[0])  # Oxygen atom connected to phosphate

    # Check if any of the phosphate oxygens are connected to the sugar carbons
    found_phosphate_on_sugar = False
    for oxygen_idx in phosphate_atoms:
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        for neighbor in oxygen_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in sugar_atoms:
                found_phosphate_on_sugar = True
                break
        if found_phosphate_on_sugar:
            break

    if not found_phosphate_on_sugar:
        return False, "Phosphate ester group not attached to monosaccharide unit"

    return True, "Contains monosaccharide unit with phosphate ester group attached to hydroxyl group"

__metadata__ = {
    'chemical_class': {
        'name': 'phospho sugar',
        'definition': 'Any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.'
    }
}