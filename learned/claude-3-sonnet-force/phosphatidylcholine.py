"""
Classifies: CHEBI:64482 phosphatidylcholine
"""
"""
Classifies: CHEBI:17719 phosphatidylcholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylcholine(smiles: str):
    """
    Determines if a molecule is a phosphatidylcholine based on its SMILES string.
    A phosphatidylcholine is a glycerophosphocholine bearing two acyl substituents at positions 1 and 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylcholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring info and identify potential glycerol backbone
    ring_info = mol.GetRingInfo()
    ring_atoms = [ring_info.AtomRings()[i] for i in range(ring_info.NumRings())]
    potential_backbone = []
    for ring in ring_atoms:
        if len(ring) == 3:
            potential_backbone.extend(ring)
    
    # Check for glycerol backbone pattern
    backbone_smarts = Chem.MolFromSmarts("[OX2]C([OX2])C([OX2])")
    backbone_match = mol.GetSubstructMatches(backbone_smarts)
    if not backbone_match:
        return False, "No glycerol backbone found"
    
    # Check for choline group
    choline_smarts = Chem.MolFromSmarts("[N+]([C])([C])([C])")
    choline_match = mol.GetSubstructMatches(choline_smarts)
    if not choline_match:
        return False, "No choline group found"
    
    # Check for phosphate group
    phosphate_smarts = Chem.MolFromSmarts("P(=O)([OX2])([OX2])[OX2]")
    phosphate_match = mol.GetSubstructMatches(phosphate_smarts)
    if not phosphate_match:
        return False, "No phosphate group found"
    
    # Check for presence of two ester groups
    ester_smarts = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_smarts)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"
    
    # Check for fatty acid chains
    chain_smarts = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_smarts)
    if len(chain_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(chain_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"
    
    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for phosphatidylcholine"
    if o_count < 6:
        return False, "Too few oxygens for phosphatidylcholine"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen (choline group)"
    
    return True, "Contains glycerol backbone with two fatty acid chains and a phosphocholine group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17719',
                          'name': 'phosphatidylcholine',
                          'definition': 'A glycerophosphocholine that is '
                                        'glycero-3-phosphocholine bearing two '
                                        'acyl substituents at positions 1 and '
                                        '2.',
                          'parents': ['CHEBI:17693', 'CHEBI:17854']},
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 147,
    'num_false_positives': 0,
    'num_true_negatives': 182445,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9866321243523316,
    'f1': 0.9932785885966384,
    'accuracy': 0.9999906878950462}