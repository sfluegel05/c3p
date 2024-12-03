"""
Classifies: CHEBI:26764 steroid hormone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_steroid_hormone(smiles: str):
    """
    Determines if a molecule is a steroid hormone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid hormone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Steroid hormones typically have a tetracyclic ring system with three 6-membered rings and one 5-membered ring
    rings = mol.GetRingInfo()
    atom_rings = rings.AtomRings()

    # Check for exactly 4 rings
    if len(atom_rings) != 4:
        return False, "Molecule does not have exactly 4 rings"

    # Check for three 6-membered rings and one 5-membered ring
    ring_sizes = [len(ring) for ring in atom_rings]
    if ring_sizes.count(6) != 3 or ring_sizes.count(5) != 1:
        return False, "Molecule does not have three 6-membered rings and one 5-membered ring"

    # Check if the rings are fused in a specific pattern
    fused_rings = [set(ring) for ring in atom_rings]
    if not all(fused_rings[i].intersection(fused_rings[i+1]) for i in range(3)):
        return False, "Rings are not fused in the typical steroid pattern"

    # Check for functional groups that are typical of hormones (e.g., hydroxyl, carbonyl, etc.)
    functional_groups = {"O": 0, "C=O": 0}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "O":
            functional_groups["O"] += 1
        if atom.GetSymbol() == "C":
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == "O" and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2:
                    functional_groups["C=O"] += 1

    if functional_groups["O"] < 1 or functional_groups["C=O"] < 1:
        return False, "Molecule does not have the necessary functional groups typical of steroid hormones"

    return True, "Molecule is classified as a steroid hormone"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26764',
                          'name': 'steroid hormone',
                          'definition': 'Any steroid that acts as hormone.',
                          'parents': ['CHEBI:35341']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '[21:06:03] SMILES Parse Error: syntax error while parsing: '
             'C\x01([C@H](C[C@@H](C/C1=C/C=C/2\\CCC[C@]3([C@]2(CC[C@@]3([C@]4(C)C[C@H](CO4)CC(O)(C)C)[H])[H])C)O)O)=C\n'
             '[21:06:03] SMILES Parse Error: Failed parsing SMILES '
             "'C\x01([C@H](C[C@@H](C/C1=C/C=C/2\\CCC[C@]3([C@]2(CC[C@@]3([C@]4(C)C[C@H](CO4)CC(O)(C)C)[H])[H])C)O)O)=C' "
             'for input: '
             "'C\x01([C@H](C[C@@H](C/C1=C/C=C/2\\CCC[C@]3([C@]2(CC[C@@]3([C@]4(C)C[C@H](CO4)CC(O)(C)C)[H])[H])C)O)O)=C'\n",
    'stdout': '',
    'num_true_positives': 3,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 26,
    'precision': 1.0,
    'recall': 0.10344827586206896,
    'f1': 0.1875,
    'accuracy': None}