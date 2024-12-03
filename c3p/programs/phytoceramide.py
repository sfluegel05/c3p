"""
Classifies: CHEBI:139051 phytoceramide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_phytoceramide(smiles: str):
    """
    Determines if a molecule is a phytoceramide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytoceramide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a fatty acid amide bond
    amide_bond = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if (atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'N') or (atom1.GetSymbol() == 'N' and atom2.GetSymbol() == 'C'):
                # Check for adjacent carbonyl group to confirm amide
                for neighbor in atom1.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom1.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        amide_bond = True
                        break
                for neighbor in atom2.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom2.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        amide_bond = True
                        break
        if amide_bond:
            break

    if not amide_bond:
        return False, "No fatty acid amide bond found"

    # Check for the presence of a phytosphingoid base
    phytosphingoid_base = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                carbons = [n for n in neighbors if n.GetSymbol() == 'C']
                if len(carbons) == 2:
                    hydroxyl_groups = 0
                    for carbon in carbons:
                        for neighbor in carbon.GetNeighbors():
                            if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                                hydroxyl_groups += 1
                    if hydroxyl_groups >= 2:
                        phytosphingoid_base = True
                        break

    if not phytosphingoid_base:
        return False, "No phytosphingoid base found"

    return True, "Molecule is a phytoceramide"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:139051',
                          'name': 'phytoceramide',
                          'definition': 'A class of phytoceramides obtained by '
                                        'formal condensation of the carboxy '
                                        'group of any fatty acid with the '
                                        'amino group of any phytosphingoid '
                                        'base.',
                          'parents': ['CHEBI:83273', 'CHEBI:84403']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 16,
    'num_false_negatives': 16,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}