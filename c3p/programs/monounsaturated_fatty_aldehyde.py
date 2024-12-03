"""
Classifies: CHEBI:61870 monounsaturated fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty aldehyde.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehyde group
    if not any(atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 2 and 
               any(neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 1 for neighbor in atom.GetNeighbors())
               for atom in mol.GetAtoms()):
        return False, "No aldehyde group found"

    # Check for exactly one double or triple bond in the aliphatic chain
    double_triple_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE)]
    
    if len(double_triple_bonds) != 1:
        return False, f"Found {len(double_triple_bonds)} double or triple bonds, expected 1"

    # Check if the chain is aliphatic
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic atoms"

    return True, "Molecule is a monounsaturated fatty aldehyde"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61870',
                          'name': 'monounsaturated fatty aldehyde',
                          'definition': 'Any fatty aldehyde with one double or '
                                        'triple bond located at any position '
                                        'in the aliphatic chain.',
                          'parents': ['CHEBI:231547']},
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
    'num_true_negatives': 10,
    'num_false_negatives': 10,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}