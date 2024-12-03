"""
Classifies: CHEBI:75885 2-pyranones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_2_pyranones(smiles: str):
    """
    Determines if a molecule is a 2-pyranone (based on the structure of 2H-pyran-2-one and its substituted derivatives).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-pyranone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all 6-membered rings with oxygen atoms
    pyranone_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if any(atom.GetSymbol() == 'O' for atom in atoms):
                pyranone_rings.append(ring)

    if not pyranone_rings:
        return False, "No 6-membered rings with oxygen found"

    # Check for the presence of a carbonyl group (C=O) attached to the ring
    for ring in pyranone_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        for atom in atoms:
            if atom.GetSymbol() == 'C':
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetIsAromatic() == False:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            return True, "2-pyranone structure detected"

    return False, "No carbonyl group (C=O) attached to the ring found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:75885',
                          'name': '2-pyranones',
                          'definition': 'A pyranone based on the structure of '
                                        '2H-pyran-2-one and its substituted '
                                        'derivatives.',
                          'parents': ['CHEBI:18946', 'CHEBI:37963']},
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
    'num_true_positives': 16,
    'num_false_positives': 15,
    'num_true_negatives': 1,
    'num_false_negatives': 0,
    'precision': 0.5161290322580645,
    'recall': 1.0,
    'f1': 0.6808510638297872,
    'accuracy': None}