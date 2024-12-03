"""
Classifies: CHEBI:51718 alpha,beta-unsaturated aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alpha_beta_unsaturated_aldehyde(smiles: str):
    """
    Determines if a molecule is an alpha,beta-unsaturated aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha,beta-unsaturated aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of aldehyde group (C=O)
    aldehyde = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        # Check for hydrogen attached to the carbon (aldehyde)
                        for neighbor2 in atom.GetNeighbors():
                            if neighbor2.GetAtomicNum() == 1:  # Hydrogen
                                aldehyde = True
                                break
                        if aldehyde:
                            break
        if aldehyde:
            break

    if not aldehyde:
        return False, "No aldehyde group found (C=O with adjacent hydrogen)"

    # Check for conjugated unsaturated C-C bond at alpha,beta position
    conjugated = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6):
                # Check for another double bond or triple bond adjacent to this double bond
                for neighbor in begin_atom.GetNeighbors():
                    if neighbor.GetIdx() != end_atom.GetIdx() and mol.GetBondBetweenAtoms(begin_atom.GetIdx(), neighbor.GetIdx()).GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
                        conjugated = True
                        break
                for neighbor in end_atom.GetNeighbors():
                    if neighbor.GetIdx() != begin_atom.GetIdx() and mol.GetBondBetweenAtoms(end_atom.GetIdx(), neighbor.GetIdx()).GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
                        conjugated = True
                        break
        if conjugated:
            break

    if not conjugated:
        return False, "No conjugated unsaturated C-C bond found at alpha,beta position"

    return True, "Molecule is an alpha,beta-unsaturated aldehyde"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51718',
                          'name': 'alpha,beta-unsaturated aldehyde',
                          'definition': 'An aldehyde of general formula '
                                        'R(1)R(2)C=CR(3)-CH=O or RC#C-CH=O in '
                                        'which the aldehydic C=O function is '
                                        'conjugated to an unsaturated C-C bond '
                                        'at the alpha,beta position.',
                          'parents': ['CHEBI:17478']},
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