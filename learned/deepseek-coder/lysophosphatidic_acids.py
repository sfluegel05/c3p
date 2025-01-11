"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: CHEBI:16480 lysophosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    A lysophosphatidic acid has a glycerol backbone with one fatty acid chain attached via an ester bond and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lysophosphatidic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with hydroxyl groups)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Look for a phosphate group (P(=O)(O)(O)-)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=[OX1])([OX2])([OX2])")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) != 1:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need exactly 1"

    # Check if the phosphate group is attached to the glycerol backbone
    phosphate_attached = False
    for match in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(match[0])
        for neighbor in phosphate_atom.GetNeighbors():
            if neighbor.GetSymbol() == "C":
                # Check if the carbon is part of the glycerol backbone
                if mol.GetSubstructMatch(glycerol_pattern):
                    phosphate_attached = True
                    break
    if not phosphate_attached:
        return False, "Phosphate group not attached to glycerol backbone"

    # Check for fatty acid chain (long carbon chain attached to ester)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - LPAs typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for lysophosphatidic acid"

    return True, "Contains glycerol backbone with one fatty acid chain and a phosphate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16480',
                          'name': 'lysophosphatidic acid',
                          'definition': 'Any monoacylglycerol phosphate obtained by hydrolytic removal of one of the two acyl groups of any phosphatidic acid or derivatives therein.',
                          'parents': ['CHEBI:16480', 'CHEBI:16480']},
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