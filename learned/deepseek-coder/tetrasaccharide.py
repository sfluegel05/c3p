"""
Classifies: CHEBI:50126 tetrasaccharide
"""
"""
Classifies: CHEBI:63505 tetrasaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrasaccharide(smiles: str):
    """
    Determines if a molecule is a tetrasaccharide based on its SMILES string.
    A tetrasaccharide is an oligosaccharide comprising four monomeric monosaccharide units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrasaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycosidic bonds (C-O-C) specific to saccharides
    glycosidic_pattern = Chem.MolFromSmarts("[C;H1,H2][O][C;H1,H2]")
    glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
    if len(glycosidic_matches) < 3:
        return False, f"Found {len(glycosidic_matches)} glycosidic bonds, need at least 3 for a tetrasaccharide"

    # Count the number of monosaccharide units by looking for characteristic patterns
    # Each monosaccharide typically has a ring structure with multiple hydroxyl groups
    monosaccharide_pattern = Chem.MolFromSmarts("[C;H1,H2][O][C;H1,H2]")
    monosaccharide_matches = mol.GetSubstructMatches(monosaccharide_pattern)
    if len(monosaccharide_matches) < 4:
        return False, f"Found {len(monosaccharide_matches)} monosaccharide units, need at least 4 for a tetrasaccharide"

    # Check molecular weight - tetrasaccharides typically have a molecular weight > 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a tetrasaccharide"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for a tetrasaccharide"
    if o_count < 8:
        return False, "Too few oxygens for a tetrasaccharide"

    return True, "Contains four monosaccharide units linked by glycosidic bonds"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63505',
                          'name': 'tetrasaccharide',
                          'definition': 'An oligosaccharide comprising four '
                                        'monomeric monosaccharide units.',
                          'parents': ['CHEBI:63505', 'CHEBI:63505']},
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