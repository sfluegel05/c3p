"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: CHEBI:37668 rotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is a rotenoid based on its SMILES string.
    A rotenoid is characterized by a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a rotenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general core rotenoid skeleton pattern
    rotenoid_pattern = Chem.MolFromSmarts("[O]1[C@@H]2[C@@H]([C@@H]3[C@@H]([C@@H]4[C@@H]([C@@H]1[C@@H]2)[C@@H]3)[C@@H]4)")
    
    # Check if the molecule contains the rotenoid skeleton
    if not mol.HasSubstructMatch(rotenoid_pattern):
        return False, "No rotenoid skeleton found"

    # Check for the presence of at least two aromatic rings (chromene part)
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)]
    if len(aromatic_rings) < 2:
        return False, "Not enough aromatic rings for chromene structure"

    # Check for the presence of at least one partially saturated ring (tetrahydro part)
    saturated_rings = [ring for ring in mol.GetRingInfo().AtomRings() if any(not mol.GetAtomWithIdx(atom).GetIsAromatic() for atom in ring)]
    if len(saturated_rings) < 1:
        return False, "No partially saturated ring found"

    # Check for the presence of oxygen atoms in the rings (chromene and tetrahydrochromene)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.IsInRing())
    if oxygen_count < 2:
        return False, "Not enough oxygen atoms in the rings"

    # Check for the presence of a carbonyl group (C=O) in the structure
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found"

    return True, "Contains a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton with appropriate substituents"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37668',
                          'name': 'rotenoid',
                          'definition': 'Members of the class of tetrahydrochromenochromene that consists of a cis-fused tetrahydrochromeno[3,4-b]chromene skeleton and its substituted derivatives. The term was originally restricted to natural products, but is now also used to describe semi-synthetic and fully synthetic compounds.'},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
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
    'accuracy': 0.9998521228585199}