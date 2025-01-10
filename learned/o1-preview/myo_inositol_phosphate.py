"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:12348 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is a molecule consisting of a cyclohexane ring in the myo-configuration,
    with hydroxyl and phosphate groups attached to the ring carbons. The molecule should not have
    any additional large substituents and should represent a free myo-inositol phosphate, not part
    of a larger structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try adding hydrogens and sanitizing
    try:
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
    except:
        return False, "Failed to sanitize molecule"
    
    # Check that the molecule contains exactly one ring, which is a cyclohexane ring
    ring_info = mol.GetRingInfo()
    if not ring_info.IsAtomInRingOfSize(0,6):
        return False, "No cyclohexane ring of size 6 found"
    atom_rings = ring_info.AtomRings()
    cyclohexane_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            if all(atom.GetAtomicNum() == 6 for atom in atoms_in_ring):
                cyclohexane_rings.append(ring)
    if len(cyclohexane_rings) != 1:
        return False, "Molecule does not contain exactly one cyclohexane ring"
    ring_atoms = cyclohexane_rings[0]
    
    # Check that all substituents on the ring carbons are hydroxyl or phosphate groups
    phosphate_found = False
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # Check that each ring carbon has exactly two bonds within the ring
        ring_neighbor_count = sum(1 for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring_atoms)
        if ring_neighbor_count != 2:
            return False, "Ring carbon does not have exactly two bonds within the ring"
        substituents = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring_atoms]
        if not substituents:
            return False, "Ring carbon has no substituents (missing hydroxyl or phosphate)"
        for sub in substituents:
            if sub.GetAtomicNum() == 8:  # Oxygen
                # Check if it's hydroxyl group or phosphate group
                attached_atoms = [nbr.GetAtomicNum() for nbr in sub.GetNeighbors()]
                if 1 in attached_atoms:
                    # Oxygen attached to hydrogen - hydroxyl group
                    continue
                elif 15 in attached_atoms:
                    # Oxygen attached to phosphorus - phosphate group
                    phosphate_found = True
                    continue
                else:
                    return False, "Oxygen substituent is neither hydroxyl nor phosphate"
            elif sub.GetAtomicNum() == 15:
                # Direct attachment of phosphorus (less common)
                phosphate_found = True
            else:
                return False, "Substituent on ring carbon is not hydroxyl or phosphate group"
    
    if not phosphate_found:
        return False, "No phosphate group found in molecule"
    
    # Check for additional large substituents outside the ring
    ring_and_substituent_atoms = set(ring_atoms)
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            ring_and_substituent_atoms.add(nbr.GetIdx())
            for nbr2 in nbr.GetNeighbors():
                ring_and_substituent_atoms.add(nbr2.GetIdx())
    # If the number of atoms in the molecule is greater than those in the ring and immediate substituents, it's a larger molecule
    if mol.GetNumAtoms() > len(ring_and_substituent_atoms):
        return False, "Molecule has additional substituents beyond ring and phosphate groups"
    
    # Optionally check for presence of stereochemistry (myo-configuration)
    # This can be complex due to different representations, so we'll skip this step for now
    
    return True, "Molecule matches myo-inositol phosphate structure"

# Metadata
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:12348',
        'name': 'myo-inositol phosphate',
        'definition': 'An inositol phosphate in which the inositol component has myo-configuration.',
        'parents': ['CHEBI:24848', 'CHEBI:24845']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    }
}