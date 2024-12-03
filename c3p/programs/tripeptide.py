"""
Classifies: CHEBI:47923 tripeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors


def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide (three amino-acid residues connected by peptide linkages).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tripeptide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for peptide bonds (C(=O)-NH)
    peptide_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'N' and 
                begin_atom.GetTotalNumHs() == 0 and end_atom.GetTotalNumHs() == 1):
                # Check if the carbon is part of a carbonyl group (C=O)
                for neighbor in begin_atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(begin_atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        peptide_bonds += 1
                        break
            elif (begin_atom.GetSymbol() == 'N' and end_atom.GetSymbol() == 'C' and 
                  begin_atom.GetTotalNumHs() == 1 and end_atom.GetTotalNumHs() == 0):
                # Check if the carbon is part of a carbonyl group (C=O)
                for neighbor in end_atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(end_atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        peptide_bonds += 1
                        break

    if peptide_bonds != 2:
        return False, f"Expected 2 peptide bonds, found {peptide_bonds}"

    # Check for three alpha amino acids
    amino_acids = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() > 0:
            # Check if the nitrogen is part of an amino group (NH2 or NH)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    # Check if the carbon is alpha (connected to a carbonyl group)
                    for n_neighbor in neighbor.GetNeighbors():
                        if n_neighbor.GetSymbol() == 'C':
                            for nn_neighbor in n_neighbor.GetNeighbors():
                                if nn_neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(n_neighbor.GetIdx(), nn_neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    amino_acids += 1
                                    break

    if amino_acids != 3:
        return False, f"Expected 3 alpha amino acids, found {amino_acids}"

    return True, "Valid tripeptide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47923',
                          'name': 'tripeptide',
                          'definition': 'Any oligopeptide that consists of '
                                        'three amino-acid residues connected '
                                        'by peptide linkages.',
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
             'position 28-29: malformed \\N character escape (<string>, line '
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