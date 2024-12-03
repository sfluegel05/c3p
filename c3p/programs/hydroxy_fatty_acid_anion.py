"""
Classifies: CHEBI:59835 hydroxy fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate anion group [C(=O)[O-]]
    carboxylate_anion = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetFormalCharge() == 0:
            neighbors = atom.GetNeighbors()
            oxygens = [n for n in neighbors if n.GetSymbol() == 'O' and n.GetFormalCharge() == -1]
            if len(oxygens) == 1:
                double_bonded_oxygens = [n for n in neighbors if n.GetSymbol() == 'O' and n.GetFormalCharge() == 0 and mol.GetBondBetweenAtoms(atom.GetIdx(), n.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE]
                if len(double_bonded_oxygens) == 1:
                    carboxylate_anion = True
                    break
    if not carboxylate_anion:
        return False, "No carboxylate anion group found"

    # Check for hydroxy group [C-OH]
    hydroxy_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetFormalCharge() == 0:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetSymbol() == 'C' and neighbors[0].GetFormalCharge() == 0:
                hydroxy_group = True
                break
    if not hydroxy_group:
        return False, "No hydroxy group found"

    return True, "Molecule is a hydroxy fatty acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:59835',
                          'name': 'hydroxy fatty acid anion',
                          'definition': 'The conjugate base of any hydroxy '
                                        'fatty acid, formed by deprotonation '
                                        'of the carboxylic acid moiety.',
                          'parents': ['CHEBI:28868', 'CHEBI:36059']},
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
    'num_true_positives': 40,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.6666666666666666,
    'recall': 1.0,
    'f1': 0.8,
    'accuracy': None}