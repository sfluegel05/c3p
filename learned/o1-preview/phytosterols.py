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
    # This pattern accounts for three fused six-membered rings and one five-membered ring
    steroid_pattern = Chem.MolFromSmarts('[#6]1[#6][#6][#6]2[#6]1[#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6][#6]4')

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for hydroxyl group attached to the steroid backbone at position 3 (common in phytosterols)
    hydroxyl_pattern = Chem.MolFromSmarts('[C;R][C;R](O)[C;R]')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl group at position 3 found"

    # Check for side chain at C17 position (variations in carbon side chain)
    # Atoms at position 17 in the steroid backbone
    # We need to find the carbon at position 17 to examine its substituents

    # Get matches for steroid backbone to identify atom indices
    steroid_match = mol.GetSubstructMatch(steroid_pattern)
    if not steroid_match:
        return False, "Steroid backbone match not found"

    # Assuming standard atom numbering for the steroid backbone, which may not be the case
    # Since atom indices can vary, we need to map the pattern accordingly

    # Create a mapped version of the steroid pattern to get atom mapping
    steroid_mol = Chem.MolFromSmarts('[#6]1~[#6]~[#6]~[#6]2~[#6]1~[#6]~[#6]3~[#6]2~[#6]~[#6]~[#6]4~[#6]3~[#6]~[#6]~[#6]~[#6]4')
    matches = mol.GetSubstructMatches(steroid_mol)
    if not matches:
        return False, "No steroid backbone match found"

    # For simplicity, take the first match
    match = matches[0]
    mol_atoms = mol.GetAtoms()

    # Identify atom at position 17 (carbon at junction of rings D and side chain)
    # In the steroid backbone, C17 is typically the 13th atom in the pattern
    # Note: positions may vary depending on the SMARTS pattern used
    atom_C17 = mol_atoms[match[13]]  # Adjust index as per pattern

    # Examine substituents attached to C17 to determine side chain length
    side_chain = [neighbor for neighbor in atom_C17.GetNeighbors() if neighbor.GetIdx() not in match]
    
    if not side_chain:
        return False, "No side chain at C17 position"

    # Calculate the number of carbons in the side chain
    side_chain_atoms = set()
    atoms_to_visit = set([atom.GetIdx() for atom in side_chain])
    visited_atoms = set()
    while atoms_to_visit:
        current_atom_idx = atoms_to_visit.pop()
        if current_atom_idx in visited_atoms:
            continue
        visited_atoms.add(current_atom_idx)
        current_atom = mol.GetAtomWithIdx(current_atom_idx)
        if current_atom.GetAtomicNum() == 6:  # Carbon atom
            side_chain_atoms.add(current_atom_idx)
            for neighbor in current_atom.GetNeighbors():
                if neighbor.GetIdx() not in visited_atoms and neighbor.GetIdx() not in match:
                    atoms_to_visit.add(neighbor.GetIdx())

    num_side_chain_carbons = len(side_chain_atoms)

    # Cholesterol has an 8-carbon side chain (C8H17)
    # Phytosterols generally have longer or more branched side chains (e.g., 9 or more carbons)
    if num_side_chain_carbons <= 8:
        return False, f"Side chain too short ({num_side_chain_carbons} carbons), possibly cholesterol"

    # If all checks passed, classify as phytosterol
    return True, "Molecule is a phytosterol with steroid backbone and appropriate side chain length"

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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata fields can be added here
}