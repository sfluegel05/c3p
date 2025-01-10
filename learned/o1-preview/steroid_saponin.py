"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:26397 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    A steroid saponin is any saponin derived from a hydroxysteroid.
    This function checks for the presence of a steroid nucleus (tetracyclic ring system),
    hydroxyl groups on the steroid nucleus, and sugar moieties attached via glycosidic bonds.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid saponin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the steroid nucleus (cyclopentanoperhydrophenanthrene)
    steroid_smarts = '[#6]1([#6])[#6][#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6]([#6])[#6]4[#6]([#6]3)[#6][#6][#6][#6]4'
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)

    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    # Check if steroid nucleus is present
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Get the atoms in the steroid nucleus
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    steroid_atoms = set()
    for match in steroid_matches:
        steroid_atoms.update(match)

    # Check for hydroxyl groups on the steroid nucleus
    hydroxyl_on_steroid = False
    for atom_idx in steroid_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check if oxygen is part of a hydroxyl group (OH)
                if neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 1:
                    hydroxyl_on_steroid = True
                    break
        if hydroxyl_on_steroid:
            break

    if not hydroxyl_on_steroid:
        return False, "No hydroxyl groups on steroid nucleus found"

    # Define SMARTS patterns for sugar moieties (pyranose and furanose rings)
    sugar_smarts_list = [
        '[OX2H][CX4H][CX4H][OX2][CX4H][CX4H]',  # Furanose ring
        '[OX2H][CX4H][CX4H][CX4H][CX4H][OX2H]'   # Pyranose ring
    ]
    sugar_patterns = [Chem.MolFromSmarts(s) for s in sugar_smarts_list]

    # Check for sugar moieties attached via glycosidic bonds
    glycosidic_bond = False
    for bond in mol.GetBonds():
        # Look for C-O single bonds
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if ((atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8) or
                (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6)):
                # Check if one atom is in the steroid nucleus and the other is in a sugar moiety
                atom1_in_steroid = atom1.GetIdx() in steroid_atoms
                atom2_in_steroid = atom2.GetIdx() in steroid_atoms

                if atom1_in_steroid != atom2_in_steroid:
                    # Identify the non-steroid atom
                    sugar_atom = atom1 if not atom1_in_steroid else atom2
                    # Check if the non-steroid atom is part of a sugar ring
                    for sugar_pattern in sugar_patterns:
                        if mol.HasSubstructMatch(sugar_pattern, useChirality=False):
                            glycosidic_bond = True
                            break
                if glycosidic_bond:
                    break

    if not glycosidic_bond:
        return False, "No sugar moieties attached via glycosidic bonds found"

    return True, "Contains steroid nucleus with hydroxyl groups and sugar moieties attached via glycosidic bonds"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:26397',
        'name': 'steroid saponin',
        'definition': 'Any saponin derived from a hydroxysteroid.',
        'parents': ['CHEBI:26023', 'CHEBI:35188']
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
    'stdout': None
}