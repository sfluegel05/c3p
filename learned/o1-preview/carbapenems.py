"""
Classifies: CHEBI:46633 carbapenems
"""
"""
Classifies: CHEBI:23066 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem is a beta-lactam antibiotic with a carbapenem skeleton,
    which is a fused beta-lactam ring (4-membered ring with nitrogen and carbonyl)
    and a five-membered ring with an exocyclic double bond,
    variously substituted at positions 3, 4, and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find rings
    ssr = Chem.GetSymmSSSR(mol)
    if len(ssr) < 2:
        return False, "Does not contain fused ring system"

    # Identify beta-lactam ring (4-membered ring with N and C=O)
    beta_lactam_ring_idx = -1
    for i, ring in enumerate(ssr):
        if len(ring) == 4:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_nums = [atom.GetAtomicNum() for atom in atoms_in_ring]
            if 7 in atom_nums:  # Nitrogen present
                # Check for carbonyl group in the ring
                has_carbonyl = False
                for bond in mol.GetBonds():
                    begin_idx = bond.GetBeginAtomIdx()
                    end_idx = bond.GetEndAtomIdx()
                    if begin_idx in ring and end_idx in ring:
                        bond_type = bond.GetBondType()
                        if bond_type == Chem.rdchem.BondType.DOUBLE:
                            begin_atom = bond.GetBeginAtom()
                            end_atom = bond.GetEndAtom()
                            if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or \
                               (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6):
                                has_carbonyl = True
                                break
                if has_carbonyl:
                    beta_lactam_ring_idx = i
                    break
    if beta_lactam_ring_idx == -1:
        return False, "Does not contain beta-lactam ring with nitrogen and carbonyl group"

    # Identify fused five-membered ring
    fused = False
    five_membered_ring_idx = -1
    for i, ring in enumerate(ssr):
        if i != beta_lactam_ring_idx and len(ring) == 5:
            # Check if rings share at least two atoms
            shared_atoms = set(ssr[beta_lactam_ring_idx]).intersection(set(ring))
            if len(shared_atoms) >= 2:
                fused = True
                five_membered_ring_idx = i
                break
    if not fused:
        return False, "Does not contain five-membered ring fused to beta-lactam ring"

    # Check for exocyclic double bond attached to the five-membered ring
    five_membered_ring_atoms = ssr[five_membered_ring_idx]
    exocyclic_double_bond_found = False
    for atom_idx in five_membered_ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in five_membered_ring_atoms:
                bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    exocyclic_double_bond_found = True
                    break
        if exocyclic_double_bond_found:
            break
    if not exocyclic_double_bond_found:
        return False, "Does not contain exocyclic double bond on five-membered ring"

    return True, "Contains carbapenem core structure"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:23066',
        'name': 'carbapenems',
        'definition': 'The class of beta-lactam antibiotics whose members have a carbapenem skeleton which is variously substituted at positions 3, 4, and 6.',
        'parents': ['CHEBI:19258', 'CHEBI:26836']
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
    'attempt': 1,
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