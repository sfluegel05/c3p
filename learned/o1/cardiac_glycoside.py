"""
Classifies: CHEBI:83970 cardiac glycoside
"""
"""
Classifies: cardiac glycoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    A cardiac glycoside has a steroid nucleus with an unsaturated lactone ring at C17
    and sugar residues attached at C3.

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

    # Define a more general steroid nucleus pattern (three 6-membered and one 5-membered rings fused)
    steroid_nucleus_smarts = '[#6&R]1-[#6&R]-[#6&R]-[#6&R]2-[#6&R]-[#6&R]-[#6&R]3-[#6&R]-[#6&R]-[#6&R]-[#6&R]-[#6&R]-[#6&R]3-[#6&R]-[#6&R]-[#6&R]-2-[#6&R]-1'
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus found"

    # Define unsaturated lactone ring (butenolide)
    lactone_smarts = 'C=1C=COC(=O)1'  # Five-membered lactone ring with unsaturation
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No unsaturated lactone ring (butenolide) found"

    # Check if lactone is attached to steroid nucleus
    # Get the matching atoms for steroid nucleus and lactone ring
    steroid_matches = mol.GetSubstructMatch(steroid_nucleus)
    lactone_matches = mol.GetSubstructMatch(lactone_pattern)

    # Check for connection between steroid nucleus and lactone ring
    connected = False
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if (begin_idx in steroid_matches and end_idx in lactone_matches) or \
           (end_idx in steroid_matches and begin_idx in lactone_matches):
            connected = True
            break
    if not connected:
        return False, "Lactone ring not attached to steroid nucleus"

    # Define sugar residue pattern (pyranose ring)
    sugar_pattern = Chem.MolFromSmarts('OC1COC(O)C(O)C1O')  # Simplified pyranose ring
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar residues found"

    # Check if sugar is attached to steroid nucleus
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    attached = False
    for sugar_match in sugar_matches:
        for bond in mol.GetBonds():
            begin_idx = bond.GetBeginAtomIdx()
            end_idx = bond.GetEndAtomIdx()
            if (begin_idx in steroid_matches and end_idx in sugar_match) or \
               (end_idx in steroid_matches and begin_idx in sugar_match):
                attached = True
                break
    if not attached:
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