"""
Classifies: CHEBI:22888 biphenyls
"""
from rdkit import Chem

def is_biphenyls(smiles: str):
    """
    Determines if a molecule is a biphenyl (benzenoid aromatic compound containing two phenyl or substituted-phenyl groups joined by a single bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biphenyl, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Find all aromatic rings
    aromatic_rings = [ring for ring in rings.AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]

    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"

    # Check for two phenyl rings connected by a single bond
    for ring1 in aromatic_rings:
        for ring2 in aromatic_rings:
            if ring1 != ring2:
                common_bonds = 0
                for atom1 in ring1:
                    for atom2 in ring2:
                        bond = mol.GetBondBetweenAtoms(atom1, atom2)
                        if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                            common_bonds += 1
                if common_bonds == 1:
                    return True, "Biphenyl structure found"

    return False, "No biphenyl structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22888',
                          'name': 'biphenyls',
                          'definition': 'Benzenoid aromatic compounds '
                                        'containing two phenyl or '
                                        'substituted-phenyl groups which are '
                                        'joined together by a single bond.',
                          'parents': [   'CHEBI:33836',
                                         'CHEBI:36820',
                                         'CHEBI:64459']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 123,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 1,
    'precision': 0.9534883720930233,
    'recall': 0.9919354838709677,
    'f1': 0.9723320158102766,
    'accuracy': None}