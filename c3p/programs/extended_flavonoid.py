"""
Classifies: CHEBI:71037 extended flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_extended_flavonoid(smiles: str):
    """
    Determines if a molecule is an extended flavonoid (flavonoid with one or more rings fused on to the phenyl substituted benzopyran framework).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an extended flavonoid, False otherwise
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

    # Find all aromatic 6-membered rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)

    if not aromatic_rings:
        return False, "No aromatic 6-membered rings found"

    # Check for flavonoid core structure (2-phenylchromen-4-one)
    flavonoid_core = Chem.MolFromSmarts('c1ccccc1-c2ccoc3ccccc23')
    if not mol.HasSubstructMatch(flavonoid_core):
        return False, "No flavonoid core structure found"

    # Check for additional fused rings
    flavonoid_core_matches = mol.GetSubstructMatches(flavonoid_core)
    for match in flavonoid_core_matches:
        core_atoms = set(match)
        fused_rings = [ring for ring in rings.AtomRings() if core_atoms.intersection(ring)]
        if len(fused_rings) > 3:  # More than three rings indicates additional fused rings
            return True, "Extended flavonoid with additional fused rings found"

    return False, "No additional fused rings found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:71037',
                          'name': 'extended flavonoid',
                          'definition': 'Any flavonoid with one or more rings '
                                        'fused on to the phenyl substituted '
                                        'benzopyran  framework.',
                          'parents': ['CHEBI:47916']},
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
    'num_true_negatives': 12,
    'num_false_negatives': 12,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}