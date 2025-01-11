"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    Glycosaminoglycans are polysaccharides containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for polysaccharide-like structure (multiple sugar units)
    # Look for patterns like -O-C-C-O- (glycosidic bonds) and -C-C-O- (sugar rings)
    glycosidic_pattern = Chem.MolFromSmarts("[OX2][CX4][CX4][OX2]")
    sugar_ring_pattern = Chem.MolFromSmarts("[CX4][CX4][OX2]")
    
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    sugar_ring_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    
    if len(glycosidic_matches) < 2 or len(sugar_ring_matches) < 2:
        return False, "Insufficient glycosidic bonds or sugar rings for a polysaccharide"

    # Check for amino sugars (aminomonosaccharides)
    # Look for patterns like -C-NH2 or -C-NH-CO- (N-acetyl groups)
    amino_sugar_pattern = Chem.MolFromSmarts("[CX4][NX3H2]")
    n_acetyl_pattern = Chem.MolFromSmarts("[CX4][NX3]([CX3](=[OX1]))[CX4]")
    
    amino_sugar_matches = mol.GetSubstructMatches(amino_sugar_pattern)
    n_acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    
    if len(amino_sugar_matches) + len(n_acetyl_matches) < 1:
        return False, "Insufficient amino sugar residues for a glycosaminoglycan"

    # Check molecular weight - glycosaminoglycans are typically large molecules
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a glycosaminoglycan"

    # Count carbons, oxygens, and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 10:
        return False, "Too few carbons for a glycosaminoglycan"
    if o_count < 5:
        return False, "Too few oxygens for a glycosaminoglycan"
    if n_count < 1:
        return False, "Too few nitrogens for a glycosaminoglycan"

    return True, "Contains polysaccharide structure with significant amino sugar residues"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18085',
                          'name': 'glycosaminoglycan',
                          'definition': 'Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.',
                          'parents': ['CHEBI:18111', 'CHEBI:47778']},
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