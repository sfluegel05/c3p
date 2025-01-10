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
    # This pattern represents three fused six-membered rings and one five-membered ring
    steroid_smarts = '[#6]1[#6][#6][#6][#6][#6]1[#6]2[#6][#6][#6][#6][#6]2[#6]3[#6][#6][#6][#6][#6]3[#6]4[#6][#6][#6][#6][#6]4'
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

    # Check for hydroxyl groups attached to the steroid nucleus
    hydroxyl_on_steroid = False
    for atom_idx in steroid_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atoms in the steroid nucleus
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                    # Check if oxygen is part of a hydroxyl group (OH)
                    if neighbor.GetDegree() == 1 and neighbor.GetTotalNumHs() == 0:
                        hydroxyl_on_steroid = True
                        break
            if hydroxyl_on_steroid:
                break

    if not hydroxyl_on_steroid:
        return False, "No hydroxyl groups on steroid nucleus found"

    # Check for sugar moieties attached via glycosidic bonds
    # Identify glycosidic bonds: C-O-C linkage where one carbon is anomeric (attached to two oxygens)
    has_sugar = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8) or \
               (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6):
                # Check for glycosidic bond
                carbon = atom1 if atom1.GetAtomicNum() == 6 else atom2
                oxygen = atom2 if atom2.GetAtomicNum() == 8 else atom1
                # Check if the carbon is anomeric (attached to two oxygens)
                num_oxygen_neighbors = sum(1 for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 8)
                if num_oxygen_neighbors >= 2:
                    # Check if one of the oxygens is connected to the steroid nucleus
                    oxygen_in_steroid = any(nbr.GetIdx() in steroid_atoms for nbr in oxygen.GetNeighbors())
                    if oxygen_in_steroid:
                        has_sugar = True
                        break

    if not has_sugar:
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
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}