"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:85726 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    A 3-oxo-fatty acyl-CoA(4-) is an acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups
    of any 3-oxo-fatty acyl-CoA, containing a long fatty acid chain with a 3-oxo group, linked to a CoA(4-) moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA(4-) moiety
    coa_4_pattern = Chem.MolFromSmarts("C(C)(COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12)NC(=O)CCNC(=O)[C@H](O)C(=O)SCCNC(=O)")
    coa_4_match = mol.HasSubstructMatch(coa_4_pattern)
    if not coa_4_match:
        return False, "Missing CoA(4-) moiety"
    
    # Check for 3-oxo group
    oxo_pattern = Chem.MolFromSmarts("CC(=O)C(=O)")
    oxo_match = mol.HasSubstructMatch(oxo_pattern)
    if not oxo_match:
        return False, "Missing 3-oxo group"
    
    # Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCCC")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"
    
    # Check fatty acid chain length
    min_chain_length = 10  # Minimum number of carbon atoms in fatty acid chain
    for match in fatty_acid_matches:
        chain_length = 0
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 6:  # Carbon
                chain_length += 1
        if chain_length >= min_chain_length:
            break
    else:
        return False, f"Fatty acid chain too short (< {min_chain_length} carbons)"
    
    # Check connectivity between 3-oxo group, fatty acid chain, and CoA(4-) moiety
    oxo_atom_idx = next((i for i, atom in enumerate(mol.GetAtoms()) if atom.GetSymbol() == 'O' and atom.GetDegree() == 1), None)
    if oxo_atom_idx is None:
        return False, "Could not find 3-oxo oxygen atom"
    oxo_neighbor = mol.GetAtomWithIdx(list(mol.GetAtomWithIdx(oxo_atom_idx).GetNeighbors())[0].GetIdx())
    if not oxo_neighbor.GetIsAromatic() and oxo_neighbor.GetDegree() == 3:
        # 3-oxo group is part of an acyl chain
        acyl_chain_start = oxo_neighbor
        acyl_chain_end = None
        for neighbor in acyl_chain_start.GetNeighbors():
            if neighbor.GetIdx() != oxo_atom_idx:
                next_atom = mol.GetAtomWithIdx(neighbor.GetIdx())
                while next_atom.GetDegree() == 2:
                    next_atom = mol.GetAtomWithIdx(list(next_atom.GetNeighbors())[0].GetIdx())
                if next_atom.GetAtomicNum() == 16:  # Sulfur (part of CoA moiety)
                    acyl_chain_end = next_atom
                    break
        if acyl_chain_end is None:
            return False, "Could not find connection between fatty acid chain and CoA(4-) moiety"
    else:
        return False, "3-oxo group not part of an acyl chain"
    
    return True, "Contains a long fatty acid chain with a 3-oxo group, linked to a CoA(4-) moiety"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:85726',
                          'name': '3-oxo-fatty acyl-CoA(4-)',
                          'definition': 'An acyl-CoA(4-) arising from '
                                        'deprotonation of the phosphate and '
                                        'diphosphate groups of any 3-oxo-fatty '
                                        'acyl-CoA.',
                          'parents': ['CHEBI:36344']},
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
    'success': False,
    'best': False,
    'error': 'Python argument types in Mol.GetAtomWithIdx(Mol, Atom) did not match C++ signature.',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': 0,
    'accuracy': None}