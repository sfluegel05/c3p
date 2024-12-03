"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for a 5-membered ring (cyclopentane)
    if not any(len(ring) == 5 for ring in rings.AtomRings()):
        return False, "No 5-membered rings found"

    # Check for a 6-membered oxygen heterocycle
    six_membered_oxygen_heterocycle = False
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                six_membered_oxygen_heterocycle = True
                break

    if not six_membered_oxygen_heterocycle:
        return False, "No 6-membered oxygen heterocycles found"

    # Check for the presence of isoprene units
    isoprene_unit = Chem.MolFromSmarts("C=C(C)C")
    if not mol.HasSubstructMatch(isoprene_unit):
        return False, "No isoprene units found"

    # Check for cleavage of the cyclopentane ring (secoiridoids)
    cyclopentane_cleaved = False
    for ring in rings.AtomRings():
        if len(ring) == 5:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetDegree() > 2 for atom in atoms):
                cyclopentane_cleaved = True
                break

    if cyclopentane_cleaved:
        return True, "Molecule is a secoiridoid monoterpenoid"
    else:
        return True, "Molecule is an iridoid monoterpenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50563',
                          'name': 'iridoid monoterpenoid',
                          'definition': 'One of a class of monoterpenoids '
                                        'biosynthesized from isoprene and '
                                        'often intermediates in the '
                                        'biosynthesis of alkaloids. Iridoids '
                                        'usually consist of a cyclopentane '
                                        'ring fused to a six-membered oxygen '
                                        'heterocycle; cleavage of a bond in '
                                        'the cyclopentane ring gives rise to '
                                        'the subclass known as secoiridoids.',
                          'parents': ['CHEBI:25409']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[00:58:54] SMILES Parse Error: syntax error while parsing: '
             'C\x01(C2(CCC(/C1=C\\C=3C=CC(=CC3)C)([H])C2(C)C)C)=O\n'
             '[00:58:54] SMILES Parse Error: Failed parsing SMILES '
             "'C\x01(C2(CCC(/C1=C\\C=3C=CC(=CC3)C)([H])C2(C)C)C)=O' for input: "
             "'C\x01(C2(CCC(/C1=C\\C=3C=CC(=CC3)C)([H])C2(C)C)C)=O'\n",
    'stdout': '',
    'num_true_positives': 6,
    'num_false_positives': 0,
    'num_true_negatives': 13,
    'num_false_negatives': 7,
    'precision': 1.0,
    'recall': 0.46153846153846156,
    'f1': 0.631578947368421,
    'accuracy': None}