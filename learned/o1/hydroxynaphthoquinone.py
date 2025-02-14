"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: CHEBI:51793 hydroxynaphthoquinone
"""
from rdkit import Chem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is any naphthoquinone where the naphthoquinone moiety is substituted by at least one hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Merge rings into fused ring systems
    fused_ring_systems = []
    visited_atoms = set()
    for ring in atom_rings:
        ring_set = set(ring)
        if not visited_atoms.intersection(ring_set):
            # Start a new fused ring system
            fused_ring_systems.append(ring_set)
            visited_atoms.update(ring_set)
        else:
            # Add to existing fused ring system
            for frs in fused_ring_systems:
                if frs.intersection(ring_set):
                    frs.update(ring_set)
                    visited_atoms.update(ring_set)
                    break

    # Check each fused ring system for naphthoquinone characteristics
    for ring_system in fused_ring_systems:
        # Check if ring system is aromatic and has at least 8 atoms (allowing for substitutions)
        if len(ring_system) >= 8:
            is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_system)
            if not is_aromatic:
                continue

            # Count ketone groups (C=O) attached to ring carbons
            ketone_count = 0
            hydroxy_count = 0

            for idx in ring_system:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    # Check for ketone group (C=O)
                    for bond in atom.GetBonds():
                        nbr = bond.GetOtherAtom(atom)
                        if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            ketone_count += 1
                            break
                    # Check for hydroxy group (C-OH)
                    for nbr in atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE and nbr.GetDegree() == 1:
                            hydroxy_count += 1
                            break

            if ketone_count >= 2 and hydroxy_count >= 1:
                return True, "Contains naphthoquinone moiety substituted with at least one hydroxy group"

    return False, "No hydroxylated naphthoquinone moiety found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51793',
        'name': 'hydroxynaphthoquinone',
        'definition': 'Any naphthoquinone in which the naphthoquinone moiety is substituted by at least one hydroxy group.',
        'parents': ['CHEBI:36124']  # CHEBI:36124 is naphthoquinone
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
    },
    'message': None,
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}