"""
Classifies: CHEBI:231547 unsaturated fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty aldehyde (fatty aldehyde containing at least one C=C or C#C bond).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of an aldehyde group (C=O at the end of the chain)
    aldehyde = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 2:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 1:
                    aldehyde = True
                    break
        if aldehyde:
            break

    if not aldehyde:
        return False, "No aldehyde group found"

    # Check for the presence of at least one C=C or C#C bond
    unsaturated = False
    for bond in mol.GetBonds():
        if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C':
                unsaturated = True
                break

    if not unsaturated:
        return False, "No C=C or C#C bond found"

    return True, "Molecule is an unsaturated fatty aldehyde"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:231547',
                          'name': 'unsaturated fatty aldehyde',
                          'definition': 'Any fatty aldehyde containing at '
                                        'least one C=C or C#C bond.',
                          'parents': ['CHEBI:35746']},
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
    'num_true_negatives': 19,
    'num_false_negatives': 19,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}