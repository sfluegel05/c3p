"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: nucleobase analogue
"""

from rdkit import Chem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.
    They typically contain pyrimidine-like or purine-like ring systems.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    has_pyrimidine_ring = False
    has_purine_ring = False

    # Get the ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    for ring in atom_rings:
        ring_size = len(ring)
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Check if the ring is aromatic
        if all(atom.GetIsAromatic() for atom in ring_atoms):
            # For pyrimidine-like rings
            if ring_size == 6:
                nitrogen_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
                if nitrogen_count >= 2:
                    has_pyrimidine_ring = True
            # For purine-like rings (fused rings)
            elif ring_size == 5 or ring_size == 6:
                # Check for fused rings
                bonds = [mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%ring_size]) for i in range(ring_size)]
                bond_rings = [ri.NumBondRings(bond.GetIdx()) for bond in bonds]
                if max(bond_rings) > 1:
                    # Ring is part of a fused system
                    fused_atoms = set(ring)
                    for other_ring in atom_rings:
                        if other_ring != ring and len(set(other_ring).intersection(fused_atoms)) > 0:
                            # Found a fused ring system
                            combined_ring = set(ring).union(other_ring)
                            combined_ring_atoms = [mol.GetAtomWithIdx(idx) for idx in combined_ring]
                            nitrogen_count = sum(1 for atom in combined_ring_atoms if atom.GetAtomicNum() == 7)
                            if len(combined_ring) >= 9 and nitrogen_count >= 4:
                                has_purine_ring = True
                                break
    if has_purine_ring:
        return True, "Contains purine-like fused ring system"
    elif has_pyrimidine_ring:
        return True, "Contains pyrimidine-like ring system"
    else:
        return False, "Does not contain purine or pyrimidine-like ring systems"
        

__metadata__ = {   'chemical_class': {   'id': None,
                              'name': 'nucleobase analogue',
                              'definition': 'A molecule that can substitute for a normal nucleobase in nucleic acids.',
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