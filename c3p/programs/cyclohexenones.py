"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone (a six-membered alicyclic ketone having one double bond in the ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    six_membered_rings = [ring for ring in rings.AtomRings() if len(ring) == 6]
    if not six_membered_rings:
        return False, "No 6-membered rings found"

    for ring in six_membered_rings:
        ring_atoms = set(ring)
        
        # Check for a ketone (C=O) in the ring
        ketone_in_ring = False
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom_idx, neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        if atom_idx in ring_atoms:
                            ketone_in_ring = True
                            break
            if ketone_in_ring:
                break
        if not ketone_in_ring:
            continue

        # Check for exactly one double bond in the ring
        double_bond_count = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom_idx = bond.GetBeginAtomIdx()
                end_atom_idx = bond.GetEndAtomIdx()
                if begin_atom_idx in ring and end_atom_idx in ring:
                    double_bond_count += 1

        if double_bond_count == 1:
            return True, "Molecule is a cyclohexenone"

    return False, "No suitable 6-membered ring with one ketone and one double bond found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:48953',
                          'name': 'cyclohexenones',
                          'definition': 'Any six-membered alicyclic ketone '
                                        'having one double bond in the ring.',
                          'parents': ['CHEBI:36132']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 52,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 1,
    'precision': 0.9629629629629629,
    'recall': 0.9811320754716981,
    'f1': 0.9719626168224299,
    'accuracy': None}