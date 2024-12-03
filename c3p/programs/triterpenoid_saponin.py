"""
Classifies: CHEBI:61778 triterpenoid saponin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_triterpenoid_saponin(smiles: str):
    """
    Determines if a molecule is a triterpenoid saponin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid saponin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycosidic bonds (at least one ether linkage to a sugar moiety)
    glycosidic_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and bond.GetBeginAtom().GetSymbol() == 'O' and bond.GetEndAtom().GetSymbol() == 'C':
            end_atom = bond.GetEndAtom()
            if any(neigh.GetSymbol() in ['O', 'N'] for neigh in end_atom.GetNeighbors()):
                glycosidic_bond = True
                break

    if not glycosidic_bond:
        return False, "No glycosidic bond found"

    # Check for triterpenoid structure (30 carbon atoms typically in a four or five-ring structure)
    triterpenoid = False
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() >= 4:
        carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
        if carbon_count >= 30:
            triterpenoid = True

    if not triterpenoid:
        return False, "No triterpenoid structure found"

    return True, "Molecule is a triterpenoid saponin"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61778',
                          'name': 'triterpenoid saponin',
                          'definition': 'A terpene glycoside in which the '
                                        'terpene moiety is a triterpenoid.',
                          'parents': ['CHEBI:26605', 'CHEBI:61777']},
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
    'num_true_positives': 48,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.8135593220338984,
    'recall': 1.0,
    'f1': 0.897196261682243,
    'accuracy': None}