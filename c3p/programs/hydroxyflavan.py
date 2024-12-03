"""
Classifies: CHEBI:72010 hydroxyflavan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hydroxyflavan(smiles: str):
    """
    Determines if a molecule is a hydroxyflavan (a flavan with one or more hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxyflavan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one flavan structure (2-phenyl-3,4-dihydro-2H-1-benzopyran)
    flavan_found = False
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetSymbol() == 'C' for atom in atoms):
                neighbors = [atom.GetNeighbors() for atom in atoms]
                if any(len(neigh) == 3 for neigh in neighbors):  # Check for the dihydro-2H-1-benzopyran structure
                    flavan_found = True
                    break

    if not flavan_found:
        return False, "No flavan structure found"

    # Check for hydroxy groups
    hydroxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neigh.GetSymbol() == 'H' for neigh in atom.GetNeighbors()):
            hydroxy_count += 1

    if hydroxy_count > 0:
        return True, f"Hydroxyflavan with {hydroxy_count} hydroxy groups"
    else:
        return False, "No hydroxy groups found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:72010',
                          'name': 'hydroxyflavan',
                          'definition': 'A member of the class of flavans in '
                                        'which one or more ring hydrogens are '
                                        'replaced by hydroxy groups.',
                          'parents': ['CHEBI:33822', 'CHEBI:38672']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 29,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}