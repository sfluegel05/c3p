"""
Classifies: CHEBI:35213 cyclodepsipeptide
"""
from rdkit import Chem

def is_cyclodepsipeptide(smiles: str):
    """
    Determines if a molecule is a cyclodepsipeptide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclodepsipeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ring structures
    rings = mol.GetRingInfo()
    if not rings.IsInitialized():
        return False, "No ring structures found"

    # Check for the presence of both amide and ester bonds in the ring
    has_amide = False
    has_ester = False

    for bond in mol.GetBonds():
        if bond.IsInRing():
            bond_type = bond.GetBondType()
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if bond_type == Chem.rdchem.BondType.SINGLE:
                if atom1.GetAtomicNum() == 7 and atom2.GetAtomicNum() == 6 and atom2.GetDegree() == 3:
                    has_amide = True
                elif atom2.GetAtomicNum() == 7 and atom1.GetAtomicNum() == 6 and atom1.GetDegree() == 3:
                    has_amide = True
                elif atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6 and atom2.GetDegree() == 3:
                    has_ester = True
                elif atom2.GetAtomicNum() == 8 and atom1.GetAtomicNum() == 6 and atom1.GetDegree() == 3:
                    has_ester = True

    if has_amide and has_ester:
        return True, "Contains both amide and ester bonds in the ring"
    else:
        return False, "Does not contain both amide and ester bonds in the ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35213',
                          'name': 'cyclodepsipeptide',
                          'definition': 'A depsipeptide in which the amino and '
                                        'hydroxy carboxylic acid residues are '
                                        'connected in a ring.',
                          'parents': ['CHEBI:23643', 'CHEBI:25000']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "'RingInfo' object has no attribute 'IsInitialized'",
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}