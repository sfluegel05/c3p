"""
Classifies: CHEBI:33708 amino-acid residue
"""
from rdkit import Chem

def is_amino_acid_residue(smiles: str):
    """
    Determines if a molecule is an amino-acid residue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino-acid residue, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a peptide bond (C-N bond in an amide group)
    peptide_bond_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if ((begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'N') or
                (begin_atom.GetSymbol() == 'N' and end_atom.GetSymbol() == 'C')):
                # Check if the carbon is double-bonded to an oxygen (amide group)
                if any(neigh.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(neigh.GetIdx(), begin_atom.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE for neigh in begin_atom.GetNeighbors()):
                    peptide_bond_found = True
                    break
                if any(neigh.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(neigh.GetIdx(), end_atom.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE for neigh in end_atom.GetNeighbors()):
                    peptide_bond_found = True
                    break

    if peptide_bond_found:
        return True, "Amino-acid residue with peptide bond found"
    else:
        return False, "No peptide bond found, not an amino-acid residue"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33708',
                          'name': 'amino-acid residue',
                          'definition': 'When two or more amino acids combine '
                                        'to form a peptide, the elements of '
                                        'water are removed, and what remains '
                                        'of each amino acid is called an '
                                        'amino-acid residue.',
                          'parents': ['CHEBI:33247']},
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
    'num_true_positives': 25,
    'num_false_positives': 6,
    'num_true_negatives': 14,
    'num_false_negatives': 32,
    'precision': 0.8064516129032258,
    'recall': 0.43859649122807015,
    'f1': 0.5681818181818181,
    'accuracy': None}