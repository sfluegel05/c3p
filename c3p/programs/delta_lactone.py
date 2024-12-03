"""
Classifies: CHEBI:18946 delta-lactone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_delta_lactone(smiles: str):
    """
    Determines if a molecule is a delta-lactone (a lactone having a six-membered lactone ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a delta-lactone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()

    # Check for a six-membered ring containing an ester (lactone) functional group
    for ring in rings:
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            # Check if the ring contains one oxygen and five carbons
            if sum(1 for atom in atoms if atom.GetSymbol() == 'O') == 1 and sum(1 for atom in atoms if atom.GetSymbol() == 'C') == 5:
                # Check if the oxygen is part of an ester (lactone)
                for atom in atoms:
                    if atom.GetSymbol() == 'O':
                        neighbors = atom.GetNeighbors()
                        if len(neighbors) == 2:
                            if (neighbors[0].GetSymbol() == 'C' and neighbors[1].GetSymbol() == 'C' and 
                                any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE for bond in mol.GetBondsBetweenAtoms(neighbors[0].GetIdx(), atom.GetIdx())) and
                                any(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in mol.GetBondsBetweenAtoms(neighbors[1].GetIdx(), atom.GetIdx()))):
                                return True, "Contains a six-membered lactone ring"
    
    return False, "Does not contain a six-membered lactone ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18946',
                          'name': 'delta-lactone',
                          'definition': 'A lactone having a six-membered '
                                        'lactone ring.',
                          'parents': ['CHEBI:25000']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "'Mol' object has no attribute 'GetBondsBetweenAtoms'",
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}