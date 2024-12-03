"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one ring
    if not rings.NumRings():
        return False, "No rings found"

    # Check for fully conjugated cyclic dione structure
    for ring in rings.AtomRings():
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if len(ring) >= 6 and all(atom.GetIsAromatic() or atom.GetHybridization() == Chem.rdchem.HybridizationType.SP2 for atom in atoms):
            carbonyl_groups = 0
            for atom in atoms:
                if atom.GetSymbol() == 'C':
                    for neigh in atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neigh.GetIdx())
                        if neigh.GetSymbol() == 'O' and bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            carbonyl_groups += 1
                            break
            if carbonyl_groups >= 2:
                return True, "Molecule has fully conjugated cyclic dione structure"

    return False, "No fully conjugated cyclic dione structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36141',
                          'name': 'quinone',
                          'definition': 'Compounds having a fully conjugated '
                                        'cyclic dione structure, such as that '
                                        'of benzoquinones, derived from '
                                        'aromatic compounds by conversion of '
                                        'an even number of -CH= groups into '
                                        '-C(=O)- groups with any necessary '
                                        'rearrangement of double bonds '
                                        '(polycyclic and heterocyclic '
                                        'analogues are included).',
                          'parents': ['CHEBI:3992']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 133,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 11,
    'precision': 0.9925373134328358,
    'recall': 0.9236111111111112,
    'f1': 0.9568345323741008,
    'accuracy': None}