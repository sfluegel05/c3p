"""
Classifies: CHEBI:73011 germacranolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_germacranolide(smiles: str):
    """
    Determines if a molecule is a germacranolide (a sesquiterpene lactone based on germacrane skeleton).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a germacranolide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 10-membered ring (germacrane skeleton)
    if not any(len(ring) == 10 for ring in rings.AtomRings()):
        return False, "No 10-membered rings found"

    # Check for lactone functionality (cyclic ester)
    lactone_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'O' and end_atom.GetSymbol() == 'C' and
                end_atom.GetDegree() == 3 and any(nb.GetSymbol() == 'O' and nb.GetIdx() != begin_atom.GetIdx() for nb in end_atom.GetNeighbors())):
                lactone_found = True
                break
            elif (end_atom.GetSymbol() == 'O' and begin_atom.GetSymbol() == 'C' and
                  begin_atom.GetDegree() == 3 and any(nb.GetSymbol() == 'O' and nb.GetIdx() != end_atom.GetIdx() for nb in begin_atom.GetNeighbors())):
                lactone_found = True
                break

    if not lactone_found:
        return False, "No lactone functionality found"

    # Check for sesquiterpene structure (15 carbon atoms)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count < 15:
        return False, f"Carbon count less than 15: {carbon_count}"

    return True, "Molecule is a germacranolide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73011',
                          'name': 'germacranolide',
                          'definition': 'A sesquiterpene lactone based on '
                                        'germacrane skeleton.',
                          'parents': ['CHEBI:37667']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[02:34:10] SMILES Parse Error: syntax error while parsing: '
             'COC(=O)C1=C/[C@H]2O[C@H]2\\C(C)=C\\[C@H]2OC(=O)C(=C)[C@@H]2[C@H](OC(=O)[C@]2(C)O[C@@H]2C)[C@H]\x01O\n'
             '[02:34:10] SMILES Parse Error: Failed parsing SMILES '
             "'COC(=O)C1=C/[C@H]2O[C@H]2\\C(C)=C\\[C@H]2OC(=O)C(=C)[C@@H]2[C@H](OC(=O)[C@]2(C)O[C@@H]2C)[C@H]\x01O' "
             'for input: '
             "'COC(=O)C1=C/[C@H]2O[C@H]2\\C(C)=C\\[C@H]2OC(=O)C(=C)[C@@H]2[C@H](OC(=O)[C@]2(C)O[C@@H]2C)[C@H]\x01O'\n",
    'stdout': '',
    'num_true_positives': 7,
    'num_false_positives': 4,
    'num_true_negatives': 6,
    'num_false_negatives': 3,
    'precision': 0.6363636363636364,
    'recall': 0.7,
    'f1': 0.6666666666666666,
    'accuracy': None}