"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: phytosterols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general steroid backbone pattern (tetracyclic ring system)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2(CC1)CCC3(C2)CCC4(C3)CCCC4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for hydroxyl group attached to the steroid backbone
    # Find atoms in the steroid backbone
    steroid_matches = mol.GetSubstructMatch(steroid_pattern)
    steroid_atoms = set(steroid_matches)

    # Identify hydroxyl groups in the molecule
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    hydroxyl_on_steroid = False
    for match in hydroxyl_matches:
        oxygen_idx = match[0]
        # Check if the oxygen is attached to a carbon in the steroid backbone
        for neighbor in mol.GetAtomWithIdx(oxygen_idx).GetNeighbors():
            if neighbor.GetIdx() in steroid_atoms:
                hydroxyl_on_steroid = True
                break
        if hydroxyl_on_steroid:
            break

    if not hydroxyl_on_steroid:
        return False, "No hydroxyl group attached to the steroid backbone"

    # Optional: Exclude cholesterol by checking the side chain
    # Cholesterol has a specific side chain at C17
    cholesterol_side_chain = Chem.MolFromSmarts('CCC(=C)C')
    side_chain_matches = mol.GetSubstructMatches(cholesterol_side_chain)
    if side_chain_matches:
        side_chain_atoms = set()
        for match in side_chain_matches:
            side_chain_atoms.update(match)
        # Check if the side chain is attached to the steroid backbone
        attachment_found = False
        for atom_idx in side_chain_atoms:
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if neighbor.GetIdx() in steroid_atoms:
                    attachment_found = True
                    break
            if attachment_found:
                break
        if attachment_found:
            return False, "Molecule is cholesterol, not a phytosterol"

    # Check for plant origin features (e.g., additional methyl or ethyl groups)
    # This is a heuristic since we can't determine origin from structure alone
    # Look for extra carbons in the side chain compared to cholesterol
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    cholesterol_carbons = 27  # Cholesterol has 27 carbons
    if num_carbons <= cholesterol_carbons:
        return False, f"Carbon count ({num_carbons}) not greater than cholesterol"

    # If all checks passed, classify as phytosterol
    return True, "Molecule is a phytosterol with steroid backbone and variations in side chain/double bonds"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'phytosterols',
        'definition': 'Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here if needed
    },
    'message': None,
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata fields can be added here
}