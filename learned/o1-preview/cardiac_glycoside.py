"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    A cardiac glycoside has a steroid nucleus (cyclopentanoperhydrophenanthrene skeleton)
    with an unsaturated lactone ring at C17 and sugar residues attached at C3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus pattern (cyclopentanoperhydrophenanthrene skeleton)
    # This pattern represents four fused rings: A, B, C (six-membered) and D (five-membered)
    steroid_nucleus_smarts = '[#6]1[#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6][#6][#6][#6]3[#6][#6][#6]2[#6]1'
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus found"

    # Define unsaturated lactone ring (butenolide)
    lactone_smarts = 'C1=COC(=O)C1'  # Five-membered unsaturated lactone ring
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No unsaturated lactone ring (butenolide) found"

    # Check if lactone is attached to steroid nucleus
    steroid_matches = mol.GetSubstructMatches(steroid_nucleus)
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    connected = False
    for steroid_match in steroid_matches:
        steroid_atoms = set(steroid_match)
        for lactone_match in lactone_matches:
            lactone_atoms = set(lactone_match)
            # Find bonds between steroid and lactone
            for bond in mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if (begin_idx in steroid_atoms and end_idx in lactone_atoms) or \
                   (end_idx in steroid_atoms and begin_idx in lactone_atoms):
                    connected = True
                    break
            if connected:
                break
        if connected:
            break
    if not connected:
        return False, "Lactone ring not attached to steroid nucleus"

    # Define sugar residue pattern (monosaccharide units)
    sugar_smarts = '[#6]1[#6][#8][#6][#6][#8][#6]1'  # Simplified pyranose ring
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar residues found"

    # Check if sugar is attached to steroid nucleus
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    sugar_attached = False
    for steroid_match in steroid_matches:
        steroid_atoms = set(steroid_match)
        for sugar_match in sugar_matches:
            sugar_atoms = set(sugar_match)
            # Find bonds between steroid and sugar
            for bond in mol.GetBonds():
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                if (begin_idx in steroid_atoms and end_idx in sugar_atoms) or \
                   (end_idx in steroid_atoms and begin_idx in sugar_atoms):
                    sugar_attached = True
                    break
            if sugar_attached:
                break
        if sugar_attached:
            break
    if not sugar_attached:
        return False, "Sugar residues not attached to steroid nucleus"

    return True, "Contains steroid nucleus with attached lactone ring and sugar residues"

__metadata__ = {
    'chemical_class': {
        'name': 'cardiac glycoside',
        'definition': 'Steroid lactones containing sugar residues that act on the contractile force of the cardiac muscles.',
    },
    'examples': [
        # SMILES strings of example cardiac glycosides provided
        'O[C@@]12[C@]3([C@@]([C@@]4([C@](CC3)(C[C@@H](O[C@@H]5O[C@H]([C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@@H](O)[C@H]5O)C)CC4)[H])C)(CC[C@@]1([C@H](CC2)C=7COC(=O)C7)C)[H])[H]',
        # Other examples can be included here
    ],
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
    },
    'message': None,
}