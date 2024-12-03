"""
Classifies: CHEBI:25235 monomethoxybenzene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_monomethoxybenzene(smiles: str):
    """
    Determines if a molecule is a monomethoxybenzene (benzene skeleton substituted with one methoxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monomethoxybenzene, False otherwise
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

    # Check if all atoms in the aromatic ring are carbon
    for ring in aromatic_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if not all(atom.GetSymbol() == 'C' for atom in atoms):
            return False, "Ring contains non-carbon atoms"

    # Check for methoxy group (-OCH3)
    methoxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 2:
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 2 and neighbors[0].GetSymbol() == 'C' and neighbors[1].GetSymbol() == 'C':
                if neighbors[1].GetSymbol() == 'C' and neighbors[1].GetDegree() == 1:
                    methoxy_count += 1

    if methoxy_count == 1:
        return True, "Molecule is a monomethoxybenzene"
    else:
        return False, f"Molecule has {methoxy_count} methoxy groups, expected exactly 1"

# Example usage
print(is_monomethoxybenzene("COc1cc(O)cc(CCCCCCC\C=C/C\C=C/CC=C)c1"))  # Example SMILES string
print(is_monomethoxybenzene("COc1ccccc1"))  # Example SMILES string without methoxy group


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25235',
                          'name': 'monomethoxybenzene',
                          'definition': 'Compounds containing a benzene '
                                        'skeleton substituted with one methoxy '
                                        'group.',
                          'parents': ['CHEBI:51683']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'Molecule has 0 methoxy groups, expected exactly 1')\n"
              "(False, 'Molecule has 0 methoxy groups, expected exactly 1')\n",
    'num_true_positives': 16,
    'num_false_positives': 9,
    'num_true_negatives': 11,
    'num_false_negatives': 17,
    'precision': 0.64,
    'recall': 0.48484848484848486,
    'f1': 0.5517241379310344,
    'accuracy': None}