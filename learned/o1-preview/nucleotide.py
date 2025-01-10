"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: CHEBI:33504 nucleotide
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    A nucleotide is a nucleoside phosphate resulting from the condensation of 
    the 3' or 5' hydroxy group of a nucleoside with phosphoric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleotide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nucleobase SMARTS patterns (common bases and modified bases)
    nucleobase_patterns = [
        # Adenine and derivatives
        Chem.MolFromSmarts('n1c[nH]c2c1ncnc2'),
        Chem.MolFromSmarts('n1c[nH]c2c1nc[nH]c2'),  # 7-deazaadenine etc.
        # Guanine and derivatives
        Chem.MolFromSmarts('O=C1NC=NC2=C1N=C[NH]2'),
        Chem.MolFromSmarts('O=C1NC=NC2=C1N=CN2'),  # Guanine
        # Cytosine and derivatives
        Chem.MolFromSmarts('N1C=CC(=O)NC1=O'),
        Chem.MolFromSmarts('N1C=CC(=O)N=C1'),  # Cytosine
        # Uracil and derivatives
        Chem.MolFromSmarts('O=C1NC=CC(=O)N1'),
        # Thymine and derivatives
        Chem.MolFromSmarts('O=C1NC(=O)C=C1C'),
        # Other modified bases may be added here
    ]

    # Check for nucleobase
    has_nucleobase = False
    for pattern in nucleobase_patterns:
        if mol.HasSubstructMatch(pattern):
            has_nucleobase = True
            break
    if not has_nucleobase:
        return False, "No nucleobase found"

    # Define sugar SMARTS patterns (ribose, deoxyribose, and modified sugars)
    sugar_patterns = [
        # Ribose
        Chem.MolFromSmarts('C1(O)C(O)C(O)C(O)O1'),
        # Deoxyribose
        Chem.MolFromSmarts('C1(O)C(O)C(O)C(O)C1'),
        # General pentose sugars
        Chem.MolFromSmarts('OC[C@H]1O[C@H](CO)[C@@H](O)[C@H]1O'),
        # Allow for modifications on the sugar
        Chem.MolFromSmarts('C1OC([C@H])(O)[C@@H](O)C1O'),
    ]

    # Check for sugar ring
    has_sugar = False
    for pattern in sugar_patterns:
        if mol.HasSubstructMatch(pattern):
            has_sugar = True
            break
    if not has_sugar:
        return False, "No pentose sugar found"

    # Define phosphate group SMARTS patterns
    phosphate_patterns = [
        Chem.MolFromSmarts('P(=O)(O)O'),  # Monophosphate
        Chem.MolFromSmarts('P(=O)(O)OP(=O)(O)O'),  # Diphosphate
        Chem.MolFromSmarts('P(=O)(O)OP(=O)(O)OP(=O)(O)O'),  # Triphosphate
        Chem.MolFromSmarts('C1OC1P(=O)(O)O'),  # Cyclic phosphate
        Chem.MolFromSmarts('OP([O-])(=O)O'),  # Deprotonated phosphate
    ]

    # Check for phosphate group
    has_phosphate = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(pattern):
            has_phosphate = True
            break
    if not has_phosphate:
        return False, "No phosphate group found"

    # Check for connectivity between nucleobase and sugar (N-glycosidic bond)
    # Look for attachment between sugar anomeric carbon and nucleobase nitrogen
    glycosidic_bond_found = False
    anomeric_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() >= 3]
    for atom in anomeric_carbons:
        bonded_atoms = atom.GetNeighbors()
        for neighbor in bonded_atoms:
            if neighbor.GetAtomicNum() == 7:
                # Check if this nitrogen is part of a nucleobase ring
                if neighbor.IsInRing():
                    glycosidic_bond_found = True
                    break
        if glycosidic_bond_found:
            break
    if not glycosidic_bond_found:
        return False, "No N-glycosidic bond between sugar and nucleobase"

    # Check for connectivity between sugar and phosphate
    phosphate_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    sugar_oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.IsInRingSize(5)]
    phosphate_linkage_found = False
    for p_atom in phosphate_atoms:
        for neighbor in p_atom.GetNeighbors():
            if neighbor in sugar_oxygen_atoms:
                phosphate_linkage_found = True
                break
        if phosphate_linkage_found:
            break
    if not phosphate_linkage_found:
        return False, "No phosphate group attached to sugar"

    return True, "Molecule is a nucleotide with nucleobase, pentose sugar, and phosphate group"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33504',
                              'name': 'nucleotide',
                              'definition': 'A nucleoside phosphate resulting from the condensation of the 3\' or 5\' hydroxy group of a nucleoside with phosphoric acid.',
                              'parents': ['CHEBI:26161']},
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
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None}