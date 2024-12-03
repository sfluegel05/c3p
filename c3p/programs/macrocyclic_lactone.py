"""
Classifies: CHEBI:63944 macrocyclic lactone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_macrocyclic_lactone(smiles: str):
    """
    Determines if a molecule is a macrocyclic lactone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrocyclic lactone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one ring
    if not rings.AtomRings():
        return False, "No rings found"

    # Check for lactone (cyclic ester) group
    lactone_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetBeginAtom().GetSymbol() == 'C' and bond.GetEndAtom().GetSymbol() == 'O':
            carbon = bond.GetBeginAtom()
            oxygen = bond.GetEndAtom()
            for neighbor in carbon.GetNeighbors():
                if neighbor.GetIdx() != oxygen.GetIdx() and neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                    lactone_found = True
                    break
            if lactone_found:
                break

    if not lactone_found:
        return False, "No lactone group found"

    # Check if the lactone is part of a macromolecule (macrocycle)
    for ring in rings.AtomRings():
        if len(ring) >= 12:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in ring_atoms):
                return True, "Macrocyclic lactone found"

    return False, "No macrocyclic lactone found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63944',
                          'name': 'macrocyclic lactone',
                          'definition': 'Any lactone in which the cyclic '
                                        'carboxylic ester group forms a part '
                                        'of a cyclic macromolecule.',
                          'parents': ['CHEBI:25000', 'CHEBI:51026']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[14:07:17] SMILES Parse Error: syntax error while parsing: '
             '[Na+].[H][C@]12C[C@]([H])(OC(=O)[C@]3([H])O[B-]45O[C@@]([H])(C(=O)O[C@@]6([H])C[C@]([H])(O[C@@H]6C)\\C=C\\C[C@@H](O)C(C)(C)[C@]6([H])CC[C@@H](C)[C@]3(O6)O4)[C@@]3(O[C@@]([H])(CC[C@H]3C)C(C)(C)[C@H](O)C\\C=C\x01)O5)[C@@H](C)O2\n'
             '[14:07:17] SMILES Parse Error: Failed parsing SMILES '
             "'[Na+].[H][C@]12C[C@]([H])(OC(=O)[C@]3([H])O[B-]45O[C@@]([H])(C(=O)O[C@@]6([H])C[C@]([H])(O[C@@H]6C)\\C=C\\C[C@@H](O)C(C)(C)[C@]6([H])CC[C@@H](C)[C@]3(O6)O4)[C@@]3(O[C@@]([H])(CC[C@H]3C)C(C)(C)[C@H](O)C\\C=C\x01)O5)[C@@H](C)O2' "
             'for input: '
             "'[Na+].[H][C@]12C[C@]([H])(OC(=O)[C@]3([H])O[B-]45O[C@@]([H])(C(=O)O[C@@]6([H])C[C@]([H])(O[C@@H]6C)\\C=C\\C[C@@H](O)C(C)(C)[C@]6([H])CC[C@@H](C)[C@]3(O6)O4)[C@@]3(O[C@@]([H])(CC[C@H]3C)C(C)(C)[C@H](O)C\\C=C\x01)O5)[C@@H](C)O2'\n"
             '[14:07:17] SMILES Parse Error: syntax error while parsing: '
             'O=C1OCC2OC(C=3NC(CC=C(C4C(OC(C=C1CC(C)C)=O)C(C(=O)O4)(CO)C)C)=CC3)=NC2/C=C\x05/C(=C(/CCCOC)\\C)/COC5CC(C)C\n'
             '[14:07:17] SMILES Parse Error: Failed parsing SMILES '
             "'O=C1OCC2OC(C=3NC(CC=C(C4C(OC(C=C1CC(C)C)=O)C(C(=O)O4)(CO)C)C)=CC3)=NC2/C=C\x05/C(=C(/CCCOC)\\C)/COC5CC(C)C' "
             'for input: '
             "'O=C1OCC2OC(C=3NC(CC=C(C4C(OC(C=C1CC(C)C)=O)C(C(=O)O4)(CO)C)C)=CC3)=NC2/C=C\x05/C(=C(/CCCOC)\\C)/COC5CC(C)C'\n",
    'stdout': '',
    'num_true_positives': 93,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 70,
    'precision': 0.9207920792079208,
    'recall': 0.5705521472392638,
    'f1': 0.7045454545454546,
    'accuracy': None}