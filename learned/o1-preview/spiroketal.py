"""
Classifies: CHEBI:72600 spiroketal
"""
"""
Classifies: CHEBI:35567 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Identify spiro atoms (atoms that are shared by exactly two rings)
    spiro_atoms = []
    for atom_idx in range(mol.GetNumAtoms()):
        num_rings = sum([atom_idx in ring for ring in atom_rings])
        if num_rings == 2:
            spiro_atoms.append(atom_idx)

    if not spiro_atoms:
        return False, "No spiro atoms found"

    # Check each spiro atom to see if it is a ketal carbon
    for spiro_idx in spiro_atoms:
        atom = mol.GetAtomWithIdx(spiro_idx)
        if atom.GetAtomicNum() != 6:
            continue  # Spiro atom is not carbon
        if atom.GetDegree() != 4:
            continue  # Spiro carbon does not have four bonds

        # Count number of oxygen and carbon neighbors
        oxy_neighbors = 0
        carbon_neighbors = 0
        for neighbor in atom.GetNeighbors():
            atomic_num = neighbor.GetAtomicNum()
            if atomic_num == 8:
                oxy_neighbors += 1
            elif atomic_num == 6:
                carbon_neighbors += 1
            else:
                break  # Neighbor is not carbon or oxygen
        else:
            if oxy_neighbors == 2 and carbon_neighbors == 2:
                return True, f"Spiroketal detected at atom index {spiro_idx}"

    return False, "No spiroketal center found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35567',
                              'name': 'spiroketal',
                              'definition': 'A cyclic ketal in which the ketal '
                                            'carbon is the only common atom of '
                                            'two rings.',
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
    'stdout': None
}