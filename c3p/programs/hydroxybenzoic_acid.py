"""
Classifies: CHEBI:24676 hydroxybenzoic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_hydroxybenzoic_acid(smiles: str):
    """
    Determines if a molecule is a hydroxybenzoic acid (benzoic acid carrying one or more phenolic hydroxy groups on the benzene ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxybenzoic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered aromatic ring
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check if one of the aromatic rings is a benzoic acid
    benzoic_acid_group = Chem.MolFromSmarts('c1ccc(cc1)C(=O)O')
    if not mol.HasSubstructMatch(benzoic_acid_group):
        return False, "No benzoic acid group found"

    # Check for phenolic hydroxy groups on the benzene ring
    phenolic_hydroxy_group = Chem.MolFromSmarts('c1cc(O)ccc1')
    if not mol.HasSubstructMatch(phenolic_hydroxy_group):
        return False, "No phenolic hydroxy groups found on the benzene ring"

    return True, "Hydroxybenzoic acid detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24676',
                          'name': 'hydroxybenzoic acid',
                          'definition': 'Any benzoic acid carrying one or more '
                                        'phenolic hydroxy groups on the '
                                        'benzene ring.',
                          'parents': ['CHEBI:22723', 'CHEBI:33853']},
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
    'num_true_positives': 41,
    'num_false_positives': 5,
    'num_true_negatives': 15,
    'num_false_negatives': 1,
    'precision': 0.8913043478260869,
    'recall': 0.9761904761904762,
    'f1': 0.9318181818181818,
    'accuracy': None}