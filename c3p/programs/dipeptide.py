"""
Classifies: CHEBI:46761 dipeptide
"""
from rdkit import Chem

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide (contains two amino-acid residues connected by peptide linkages).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the list of peptide bonds (C-N bonds in amide groups)
    peptide_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 7) or (begin_atom.GetAtomicNum() == 7 and end_atom.GetAtomicNum() == 6):
                # Check if the bond is part of an amide group
                if (begin_atom.GetAtomicNum() == 6 and begin_atom.GetTotalDegree() == 3 and any(neighbor.GetAtomicNum() == 8 for neighbor in begin_atom.GetNeighbors())) or \
                   (end_atom.GetAtomicNum() == 6 and end_atom.GetTotalDegree() == 3 and any(neighbor.GetAtomicNum() == 8 for neighbor in end_atom.GetNeighbors())):
                    peptide_bonds.append(bond)

    if len(peptide_bonds) != 1:
        return False, "The molecule does not contain exactly one peptide bond"

    # Check if there are exactly two amino acid residues
    amino_acid_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and any(neighbor.GetAtomicNum() == 6 and neighbor.GetTotalDegree() == 3 and any(n.GetAtomicNum() == 8 for n in neighbor.GetNeighbors()) for neighbor in atom.GetNeighbors()):
            amino_acid_count += 1

    if amino_acid_count != 2:
        return False, "The molecule does not contain exactly two amino acid residues"

    return True, "The molecule is a dipeptide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46761',
                          'name': 'dipeptide',
                          'definition': 'Any molecule that contains two '
                                        'amino-acid residues connected by '
                                        'peptide linkages.',
                          'parents': ['CHEBI:25676']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 30-31: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}