"""
Classifies: CHEBI:26347 prostanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_prostanoid(smiles: str):
    """
    Determines if a molecule is a prostanoid (natural prostaglandins and prostaglandin-like compounds including prostacyclins and thromboxanes).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostanoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Prostanoids typically have a cyclopentane ring with two side chains
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 5 for ring in ring_info.AtomRings()):
        return False, "No 5-membered rings found"

    # Check for specific functional groups and structural features
    has_carboxylic_acid = any(atom.GetSymbol() == 'C' and atom.GetTotalNumHs() == 0 and 
                              any(neighbor.GetSymbol() == 'O' and neighbor.GetTotalNumHs() == 1 for neighbor in atom.GetNeighbors())
                              for atom in mol.GetAtoms())
    has_hydroxy_groups = any(atom.GetSymbol() == 'O' and atom.GetTotalNumHs() == 1 for atom in mol.GetAtoms())
    has_double_bonds = any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in mol.GetBonds())

    # Prostanoids often have a carboxylic acid group at the end of a side chain
    pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(pattern):
        return False, "No carboxylic acid group found"

    if not has_hydroxy_groups:
        return False, "No hydroxy groups found"

    if not has_double_bonds:
        return False, "No double bonds found"

    return True, "Molecule is classified as a prostanoid"

# Example usage:
smiles = "O[C@H]1[C@@H]([C@H](C(=O)C1)C/C=C\\CCCC(O)=O)/C=C/[C@@H](O)CCCCC"
is_prostanoid(smiles)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26347',
                          'name': 'prostanoid',
                          'definition': 'The family of natural prostaglandins '
                                        'and prostaglandin-like compounds '
                                        'including prostacyclins and '
                                        'thromboxanes.',
                          'parents': ['CHEBI:23899']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 51,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 10,
    'precision': 0.9622641509433962,
    'recall': 0.8360655737704918,
    'f1': 0.8947368421052632,
    'accuracy': None}