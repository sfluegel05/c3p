"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
"""

from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify quinone ring
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    quinone_ring_found = False
    quinone_ring_atoms = None

    for ring in atom_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        num_carbonyls = 0
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 6:  # Carbon
                for neighbor in atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor_idx)
                    if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        num_carbonyls +=1
                        break
        if num_carbonyls >=2:
            # Found a quinone ring
            quinone_ring_found = True
            quinone_ring_atoms = set(ring)
            break  # Stop at first quinone ring found

    if not quinone_ring_found:
        return False, "No quinone ring found"

    # Identify side chains attached to the quinone ring
    side_chains = []
    for atom_idx in quinone_ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in quinone_ring_atoms:
                # Found an atom outside the ring
                # Traverse the side chain starting from neighbor
                side_chain_atoms = set()
                atoms_to_visit = [neighbor_idx]
                while atoms_to_visit:
                    current_idx = atoms_to_visit.pop()
                    if current_idx not in side_chain_atoms and current_idx not in quinone_ring_atoms:
                        side_chain_atoms.add(current_idx)
                        current_atom = mol.GetAtomWithIdx(current_idx)
                        for nbr in current_atom.GetNeighbors():
                            nbr_idx = nbr.GetIdx()
                            if nbr_idx not in side_chain_atoms and nbr_idx not in quinone_ring_atoms:
                                atoms_to_visit.append(nbr_idx)
                side_chains.append(side_chain_atoms)

    # Analyze side chains
    prenyl_side_chain_found = False
    for side_chain_atoms in side_chains:
        num_carbons = 0
        num_double_bonds = 0
        num_methyl_branches = 0
        visited_bonds = set()
        for idx in side_chain_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                num_carbons +=1
                neighbor_indices = [a.GetIdx() for a in atom.GetNeighbors() if a.GetIdx() in side_chain_atoms]
                if len(neighbor_indices) > 2:
                    # Branch point
                    num_methyl_branches +=1
            for bond in atom.GetBonds():
                bond_idx = bond.GetIdx()
                if bond_idx in visited_bonds:
                    continue
                visited_bonds.add(bond_idx)
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if begin_idx in side_chain_atoms and end_idx in side_chain_atoms:
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        num_double_bonds +=1
        if num_carbons >=10 and num_double_bonds >=2 and num_methyl_branches >=2:
            prenyl_side_chain_found = True
            break

    if not prenyl_side_chain_found:
        return False, "No polyprenyl-derived side chain attached to quinone ring"

    return True, "Contains quinone ring with polyprenyl-derived side chain"

__metadata__ = {
    'chemical_class': {
        'id': '',
        'name': 'prenylquinone',
        'definition': 'A quinone substituted by a polyprenyl-derived side-chain. Prenylquinones occur in all living cells. Due to their amphiphilic character, they are mainly located in biological membranes where they function as electron and proton carriers in the photosynthetic and respiratory electron transport chains. Some prenylquinones also perform more specialised roles such as antioxidants and enzyme cofactors. Prenylquinones are classified according to ring structure: the main classes are menaquinones, phylloquinones, ubiquinones and plastoquinones.',
        'parents': []
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
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}