"""
Classifies: CHEBI:22723 benzoic acids
"""
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_benzoic_acids(smiles: str):
    """
    Determines if a molecule is a benzoic acid (aromatic carboxylic acid with at least one carboxy group attached to a benzene ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a benzoic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()

    # Check for at least one benzene ring
    benzene_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms) and all(atom.GetSymbol() == 'C' for atom in atoms):
                benzene_rings.append(ring)

    if not benzene_rings:
        return False, "No benzene rings found"

    # Check for carboxy group attached to benzene ring
    carboxy_groups = []
    for ring in benzene_rings:
        ring_atoms = set(ring)
        for atom_idx in ring_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 3:
                    # Check if this carbon is part of a carboxy group (C(=O)O)
                    carboxy = False
                    for n in neighbor.GetNeighbors():
                        if n.GetSymbol() == 'O' and n.GetDegree() == 1:
                            for nn in n.GetNeighbors():
                                if nn.GetSymbol() == 'C' and nn.GetIdx() == neighbor.GetIdx():
                                    carboxy = True
                                    break
                    if carboxy:
                        carboxy_groups.append((atom_idx, neighbor.GetIdx()))

    if not carboxy_groups:
        return False, "No carboxy groups attached to benzene rings found"

    return True, "Benzoic acid detected"

# Example usage
smiles = "O=C(O)C1=CC=CC=C1"  # Benzoic acid
print(is_benzoic_acids(smiles))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22723',
                          'name': 'benzoic acids',
                          'definition': 'Any aromatic carboxylic acid that '
                                        'consists of benzene in which at least '
                                        'a single hydrogen has been '
                                        'substituted by a carboxy group.',
                          'parents': ['CHEBI:22712', 'CHEBI:33859']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 16-17: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}