"""
Classifies: CHEBI:61655 steroid saponin
"""
"""
Classifies: CHEBI:26397 steroid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Kekulize the molecule to handle aromaticity
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except:
        pass

    # Define a more general steroid nucleus SMARTS pattern
    # This pattern allows for variations in saturation and ring junctions
    steroid_pattern = Chem.MolFromSmarts(
        "[#6]1([#6])[#6][#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6][#6]4[#6]([#6]3)[#6][#6][#6][#6]4"
    )

    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    # Check if steroid nucleus is present
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Find all matches to the steroid nucleus
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    if len(steroid_matches) == 0:
        return False, "No steroid nucleus found"

    # Get the atoms in the steroid nucleus
    # For simplicity, consider the first match
    steroid_atoms = set(steroid_matches[0])

    # Check for hydroxyl groups (or other oxygen functionalities) on the steroid nucleus
    hydroxyl_on_steroid = False
    for atom_idx in steroid_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Check for any oxygen atoms bonded to the steroid atoms
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                hydroxyl_on_steroid = True
                break
        if hydroxyl_on_steroid:
            break

    if not hydroxyl_on_steroid:
        return False, "No hydroxyl groups on steroid nucleus found"

    # Define a sugar ring pattern (5 or 6 membered ring with oxygen and hydroxyl groups)
    sugar_pattern = Chem.MolFromSmarts(
        "[#6&R1]-O-[#6&R1]-[#6&R1]-[#6&R1]-[#6&R1]-O"  # Example for pyranose ring
    )
    if sugar_pattern is None:
        return False, "Invalid sugar SMARTS pattern"

    # Check for sugar moieties attached via glycosidic bonds
    glycosidic_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Check if bond is between oxygen and carbon
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
               (atom2.GetAtomicNum() == 8 and atom1.GetAtomicNum() == 6):
                # Check if one atom is part of the steroid nucleus and the other is part of a sugar ring
                atom1_in_steroid = atom1.GetIdx() in steroid_atoms
                atom2_in_steroid = atom2.GetIdx() in steroid_atoms

                # Identify the sugar ring containing the atom
                if not atom1_in_steroid:
                    is_sugar = atom1.IsInRingSize(5) or atom1.IsInRingSize(6)
                else:
                    is_sugar = atom2.IsInRingSize(5) or atom2.IsInRingSize(6)

                if (atom1_in_steroid != atom2_in_steroid) and is_sugar:
                    glycosidic_bond = True
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
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None
}