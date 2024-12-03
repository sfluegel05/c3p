"""
Classifies: CHEBI:24613 homodetic cyclic peptide
"""
from rdkit import Chem

def is_homodetic_cyclic_peptide(smiles: str):
    """
    Determines if a molecule is a homodetic cyclic peptide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a homodetic cyclic peptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is cyclic
    if not mol.GetRingInfo().AtomRings():
        return False, "Molecule is not cyclic"

    # Check for peptide bonds (amide linkages)
    peptide_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'N') or (atom1.GetSymbol() == 'N' and atom2.GetSymbol() == 'C'):
                # Check if the carbon is part of a carbonyl group (C=O)
                if any(nb.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom1.GetIdx(), nb.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE for nb in atom1.GetNeighbors()) or \
                   any(nb.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom2.GetIdx(), nb.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE for nb in atom2.GetNeighbors()):
                    peptide_bond_count += 1

    if peptide_bond_count == 0:
        return False, "No peptide bonds found"

    # Check if all rings consist solely of amino-acid residues in peptide linkages
    for ring in mol.GetRingInfo().AtomRings():
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() not in ['C', 'N', 'O', 'H']:
                return False, "Ring contains non-amino acid atoms"

    return True, "Molecule is a homodetic cyclic peptide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24613',
                          'name': 'homodetic cyclic peptide',
                          'definition': 'A homodetic cyclic peptide is a '
                                        'cyclic peptide in which the ring '
                                        'consists solely of amino-acid '
                                        'residues in peptide linkages.',
                          'parents': ['CHEBI:23449']},
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
    'num_true_positives': 13,
    'num_false_positives': 7,
    'num_true_negatives': 9,
    'num_false_negatives': 3,
    'precision': 0.65,
    'recall': 0.8125,
    'f1': 0.7222222222222223,
    'accuracy': None}