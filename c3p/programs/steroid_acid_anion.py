"""
Classifies: CHEBI:50160 steroid acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_steroid_acid_anion(smiles: str):
    """
    Determines if a molecule is a steroid acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid structure: three 6-membered rings and one 5-membered ring
    steroid_rings = mol.GetRingInfo().AtomRings()
    if len(steroid_rings) < 4:
        return False, "Molecule does not have enough rings for a typical steroid structure"
    
    six_membered_rings = [ring for ring in steroid_rings if len(ring) == 6]
    five_membered_rings = [ring for ring in steroid_rings if len(ring) == 5]

    if len(six_membered_rings) < 3 or len(five_membered_rings) < 1:
        return False, "Molecule does not have the typical steroid ring structure (three 6-membered rings and one 5-membered ring)"

    # Check for carboxylate anion group (COO-)
    carboxylate_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetFormalCharge() == 0:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                o_atoms = [neighbor for neighbor in neighbors if neighbor.GetSymbol() == 'O']
                if len(o_atoms) == 2:
                    if (o_atoms[0].GetFormalCharge() == -1 and o_atoms[1].GetFormalCharge() == 0) or (o_atoms[0].GetFormalCharge() == 0 and o_atoms[1].GetFormalCharge() == -1):
                        carboxylate_found = True
                        break

    if not carboxylate_found:
        return False, "No carboxylate anion group found"

    return True, "Molecule is a steroid acid anion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50160',
                          'name': 'steroid acid anion',
                          'definition': 'Any anion formed by loss of a proton '
                                        'from a steroid acid.',
                          'parents': ['CHEBI:29067']},
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
    'num_true_positives': 15,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}