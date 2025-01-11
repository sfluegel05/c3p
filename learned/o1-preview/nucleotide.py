"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33504 nucleotide
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate_pattern = Chem.MolFromSmarts('P(=O)(O)(O)O')
    has_phosphate = mol.HasSubstructMatch(phosphate_pattern)
    if not has_phosphate:
        return False, "No phosphate group found"

    # Check for sugar (five-membered ring with oxygen)
    sugar_pattern = Chem.MolFromSmarts('C1OC[C@H](O)[C@@H]1O')  # Ribose
    has_sugar = mol.HasSubstructMatch(sugar_pattern)
    if not has_sugar:
        # Try a more general sugar pattern allowing for modifications
        sugar_pattern = Chem.MolFromSmarts('C1OC([C@H])([C@@H])O1')  # Generalized furanose ring
        has_sugar = mol.HasSubstructMatch(sugar_pattern)
        if not has_sugar:
            return False, "No pentose sugar found"

    # Check for nucleobase (aromatic ring with at least two nitrogen atoms)
    nucleobase_found = False
    for ring in mol.GetRingInfo().AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if all(atom.GetIsAromatic() for atom in ring_atoms):
            num_nitrogens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
            if num_nitrogens >= 2:
                nucleobase_found = True
                break
    if not nucleobase_found:
        return False, "No nucleobase found"

    # Check for connection between nucleobase and sugar (N-glycosidic bond)
    has_nucleoside_linkage = False
    for bond in mol.GetBonds():
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        # Look for bond between nitrogen (from nucleobase) and carbon (from sugar)
        if (begin_atom.GetAtomicNum() == 7 and end_atom.GetAtomicNum() == 6) or \
           (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 7):
            # Check if one atom is in an aromatic ring (nucleobase) and the other in sugar ring
            begin_in_base = begin_atom.GetIsAromatic()
            end_in_base = end_atom.GetIsAromatic()
            begin_in_sugar = mol.GetRingInfo().IsAtomInRingOfSize(begin_atom.GetIdx(), 5)
            end_in_sugar = mol.GetRingInfo().IsAtomInRingOfSize(end_atom.GetIdx(), 5)
            if (begin_in_base and end_in_sugar) or (begin_in_sugar and end_in_base):
                has_nucleoside_linkage = True
                break
    if not has_nucleoside_linkage:
        return False, "No nucleoside linkage between nucleobase and sugar"

    # Check for connection between sugar and phosphate
    has_phosphate_linkage = False
    phosphate_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    sugar_atoms = []
    for match in mol.GetSubstructMatches(sugar_pattern):
        sugar_atoms.extend(match)
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if (begin_idx in phosphate_atoms and end_idx in sugar_atoms) or \
           (begin_idx in sugar_atoms and end_idx in phosphate_atoms):
            has_phosphate_linkage = True
            break
    if not has_phosphate_linkage:
        return False, "No phosphate group attached to sugar"

    return True, "Molecule is a nucleotide with nucleobase, pentose sugar, and phosphate group"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33504',
                              'name': 'nucleotide',
                              'definition': 'A nucleoside phosphate resulting from the condensation of the 3\' or 5\' hydroxy group of a nucleoside with phosphoric acid.',
                              'parents': ['CHEBI:26161']},
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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None}