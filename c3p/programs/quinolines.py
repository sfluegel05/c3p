"""
Classifies: CHEBI:26513 quinolines
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_quinolines(smiles: str):
    """
    Determines if a molecule is a quinoline (benzene ring ortho fused to carbons 2 and 3 of a pyridine ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinoline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generate the ring information
    rings = mol.GetRingInfo()

    # Find all rings
    atom_rings = rings.AtomRings()

    # Check for the presence of a benzene ring and a pyridine ring
    benzene_rings = []
    pyridine_rings = []
    for ring in atom_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if len(ring) == 6 and all(atom.GetSymbol() == 'C' for atom in atoms):
            benzene_rings.append(ring)
        elif len(ring) == 6 and sum(1 for atom in atoms if atom.GetSymbol() == 'N') == 1:
            pyridine_rings.append(ring)

    if not benzene_rings:
        return False, "No benzene rings found"
    if not pyridine_rings:
        return False, "No pyridine rings found"

    # Check for ortho fusion between benzene and pyridine rings at carbons 2 and 3 of pyridine
    for benzene_ring in benzene_rings:
        for pyridine_ring in pyridine_rings:
            common_atoms = set(benzene_ring).intersection(pyridine_ring)
            if len(common_atoms) == 2:
                pyridine_atoms = [mol.GetAtomWithIdx(i) for i in pyridine_ring]
                pyridine_indices = [atom.GetIdx() for atom in pyridine_atoms]
                # Check if the common atoms are adjacent in the pyridine ring
                common_indices = sorted(common_atoms)
                if pyridine_indices.index(common_indices[1]) - pyridine_indices.index(common_indices[0]) == 1:
                    return True, "Quinoline structure found"
    
    return False, "No ortho fusion between benzene and pyridine rings at carbons 2 and 3 of pyridine"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26513',
                          'name': 'quinolines',
                          'definition': 'A class of aromatic heterocyclic '
                                        'compounds each of which contains a '
                                        'benzene ring ortho fused to carbons 2 '
                                        'and 3 of a pyridine ring.',
                          'parents': ['CHEBI:33659', 'CHEBI:38101']},
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
    'num_true_positives': 103,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 73,
    'precision': 1.0,
    'recall': 0.5852272727272727,
    'f1': 0.7383512544802867,
    'accuracy': None}