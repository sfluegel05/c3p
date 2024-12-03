"""
Classifies: CHEBI:83273 N-acylsphingoid
"""
from rdkit import Chem

def is_N_acylsphingoid(smiles: str):
    """
    Determines if a molecule is an N-acylsphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an amide bond
    amide_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'N' and end_atom.GetSymbol() == 'C' and 
                any(neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(end_atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE for neighbor in end_atom.GetNeighbors())):
                amide_bond = True
                break
            if (end_atom.GetSymbol() == 'N' and begin_atom.GetSymbol() == 'C' and 
                any(neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(begin_atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE for neighbor in begin_atom.GetNeighbors())):
                amide_bond = True
                break

    if not amide_bond:
        return False, "No amide bond found"

    # Check for the sphingoid base (long-chain amino alcohol)
    sphingoid_base = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            carbon_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'C']
            if len(carbon_neighbors) >= 1:
                for carbon in carbon_neighbors:
                    oxygen_neighbors = [n for n in carbon.GetNeighbors() if n.GetSymbol() == 'O']
                    if len(oxygen_neighbors) >= 1:
                        # Check for long chain (at least 12 carbons)
                        chain_length = 1
                        current_atom = carbon
                        visited_atoms = set()
                        while current_atom and chain_length < 12:
                            visited_atoms.add(current_atom.GetIdx())
                            next_atoms = [n for n in current_atom.GetNeighbors() if n.GetSymbol() == 'C' and n.GetIdx() not in visited_atoms]
                            if next_atoms:
                                current_atom = next_atoms[0]
                                chain_length += 1
                            else:
                                current_atom = None
                        if chain_length >= 12:
                            sphingoid_base = True
                            break
            if sphingoid_base:
                break

    if not sphingoid_base:
        return False, "No sphingoid base found"

    return True, "Molecule is an N-acylsphingoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:83273',
                          'name': 'N-acylsphingoid',
                          'definition': 'A ceramide consisting of an undefined '
                                        'sphingoid base linked to an undefined '
                                        'fatty acid via an amide bond.',
                          'parents': ['CHEBI:17761']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 28,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 11,
    'precision': 0.9333333333333333,
    'recall': 0.717948717948718,
    'f1': 0.8115942028985509,
    'accuracy': None}