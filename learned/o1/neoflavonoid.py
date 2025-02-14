"""
Classifies: CHEBI:71971 neoflavonoid
"""
"""
Classifies: neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define benzopyran core with position 4 labeled
    benzopyran_smiles = 'O1C=CC=CC=C1'  # SMILES for 1-benzopyran
    benzopyran_mol = Chem.MolFromSmiles(benzopyran_smiles)
    if benzopyran_mol is None:
        return False, "Failed to create benzopyran core"

    # Label the atom at position 4 (opposite the oxygen in the pyran ring)
    # Atom indices in benzopyran_mol:
    # 0: O, 1: C, 2: C, 3: C, 4: C, 5: C, 6: C
    atom4_idx = 4  # Zero-based index for atom at position 4
    benzopyran_mol.GetAtomWithIdx(atom4_idx).SetAtomMapNum(4)

    # Generate the SMARTS pattern
    benzopyran_smarts = Chem.MolToSmarts(benzopyran_mol)
    benzopyran_pattern = Chem.MolFromSmarts(benzopyran_smarts)
    if benzopyran_pattern is None:
        return False, "Invalid benzopyran SMARTS pattern"

    # Search for the benzopyran core in the molecule
    matches = mol.GetSubstructMatches(benzopyran_pattern, useChirality=False)
    if not matches:
        return False, "1-benzopyran core not found"

    # For each match, check if atom at position 4 has an aryl substituent
    for match in matches:
        # Find the index of the atom with atom map number 4 in the pattern
        for pat_idx, pat_atom in enumerate(benzopyran_pattern.GetAtoms()):
            if pat_atom.GetAtomMapNum() == 4:
                pattern_atom4_idx = pat_idx
                break
        else:
            continue  # Atom map 4 not found in pattern

        # Map from pattern atom index to molecule atom index
        mol_atom_idx_4 = match[pattern_atom4_idx]
        atom4 = mol.GetAtomWithIdx(mol_atom_idx_4)

        # Get neighbors of atom4 not in the benzopyran core
        benzopyran_atom_indices = set(match)
        neighbor_atoms = [a for a in atom4.GetNeighbors() if a.GetIdx() not in benzopyran_atom_indices]

        # Check if any neighbor is part of an aromatic ring (aryl group)
        for neighbor in neighbor_atoms:
            if neighbor.GetIsAromatic():
                # Check if neighbor is in an aromatic ring
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if neighbor.GetIdx() in ring:
                        # Check if ring is aromatic
                        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                            return True, "Molecule is a neoflavonoid with aryl substituent at position 4"
            # Alternatively, check if neighbor is connected to an aromatic ring
            else:
                for nbr_neighbor in neighbor.GetNeighbors():
                    if nbr_neighbor.GetIsAromatic():
                        ring_info = mol.GetRingInfo()
                        for ring in ring_info.AtomRings():
                            if nbr_neighbor.GetIdx() in ring:
                                if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                                    return True, "Molecule is a neoflavonoid with aryl substituent at position 4"

    return False, "No aryl substituent at position 4 of the 1-benzopyran core"


__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'neoflavonoid',
                              'definition': 'Any 1-benzopyran with an aryl substituent at position 4. The term was originally restricted to natural products, but is now also used to describe semi-synthetic and fully synthetic compounds.',
                              'parents': []},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None}