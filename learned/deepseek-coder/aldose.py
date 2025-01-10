"""
Classifies: CHEBI:15693 aldose
"""
"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is a polyhydroxy aldehyde (H[CH(OH)]nC(=O)H, n >= 2) or its intramolecular hemiacetal.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aldose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aldehyde group (C=O) in open-chain or cyclic form
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        # Check for hemiacetal form (cyclic structure with oxygen in the ring)
        hemiacetal_pattern = Chem.MolFromSmarts("[OX2;R][CX4;R][CX3;R](=O)")
        if not mol.HasSubstructMatch(hemiacetal_pattern):
            return False, "No aldehyde or hemiacetal group found"

    # Count hydroxyl groups (OH)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    if hydroxyl_count < 2:
        return False, f"Found {hydroxyl_count} hydroxyl groups, need at least 2"

    # Check for carbon chain length (n >= 2)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 3:
        return False, f"Found {carbon_count} carbons, need at least 3"

    # Check for cyclic structure (optional, as aldoses can be open-chain or cyclic)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Ensure the ring contains an oxygen (hemiacetal form)
        ring_atoms = set()
        for ring in ring_info.AtomRings():
            ring_atoms.update(ring)
        oxygen_in_ring = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring_atoms)
        if not oxygen_in_ring:
            return False, "Cyclic structure without oxygen (not a hemiacetal)"

    return True, "Contains aldehyde or hemiacetal group with multiple hydroxyl groups"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15693',
        'name': 'aldose',
        'definition': 'Aldehydic parent sugars (polyhydroxy aldehydes H[CH(OH)]nC(=O)H, n >= 2) and their intramolecular hemiacetals.',
        'parents': ['CHEBI:26647', 'CHEBI:47778']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199
}