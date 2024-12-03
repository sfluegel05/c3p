"""
Classifies: CHEBI:38757 isoflavones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one 6-membered ring
    if not any(len(ring) == 6 for ring in rings.AtomRings()):
        return False, "No 6-membered rings found"

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check for the presence of a 3-aryl-1-benzopyran-4-one skeleton
    # This requires a benzopyran ring system with an aryl group at position 3 and a ketone at position 4
    benzopyran_pattern = Chem.MolFromSmarts('c1cc2oc(=O)cc(c2c1)-c3ccccc3')
    if not mol.HasSubstructMatch(benzopyran_pattern):
        return False, "No 3-aryl-1-benzopyran-4-one skeleton found"

    return True, "Molecule is an isoflavone"

# Example usage:
# smiles = "COc1ccc(cc1)-c1coc2cc(O)c(OC)cc2c1=O"
# result, reason = is_isoflavones(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38757',
                          'name': 'isoflavones',
                          'definition': 'Any isoflavonoid with a '
                                        '3-aryl-1-benzopyran-4-one '
                                        '(3-aryl-4H-chromen-4-one) skeleton '
                                        'and its substituted derivatives.',
                          'parents': ['CHEBI:3992', 'CHEBI:50753']},
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
    'num_false_negatives': 25,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}